#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
from __future__ import annotations

import re
from typing import Generator, Iterable, Iterator, Optional, Sequence

from paleomix.common.formats import FormatError
from paleomix.common.utilities import Immutable, TotallyOrdered

NodeID = Optional[int]
NodeName = Optional[str]


class NewickError(FormatError):
    pass


class NewickParseError(NewickError):
    """Exception raised if errors occur during parsing
    of Newick strings."""


class Newick(TotallyOrdered, Immutable):
    """Immutable object representing a Newick node.

    Nodes are classified as either internal nodes (have children),
    or leaf nodes (does not have children). A node MUST either have
    1 or more child-nodes, or have a name and/or a length. This is to
    ensure that nodes can be represented in an unambigious manner
    using the Newick format.

    No assumptions are made about the type of the 'name' and the 'length'
    properties when simply parsing the tree, and these are simply converted
    into strings when the Newick string is generated. However, additional
    contraints apply when unrooting/rerooting trees (see below).
    """

    name: str | None
    length: float | str | None
    children: tuple[Newick, ...]
    _hash: int
    _weight: int

    def __init__(
        self,
        name: str | None = None,
        length: str | float | None = None,
        children: Iterable[Newick] = (),
    ) -> None:
        children = tuple(children)
        name = name or None

        Immutable.__init__(
            self,
            name=name,
            length=length,
            children=children,
            _hash=hash((name, length, children)),
        )

        if not (self.children or self.name is not None or self.length is not None):
            raise NewickError("Leaf nodes MUST have either a name or a length")

        weight = 0
        for child in self.children:
            if not isinstance(child, Newick):
                raise TypeError("Child nodes must be Newick nodes")
            weight += 1
        object.__setattr__(self, "_weight", weight)

    @property
    def is_leaf(self) -> bool:
        """Returns true if the node is a leaf (has no children)."""
        return not self.children

    def get_leaf_nodes(self) -> Iterator[Newick]:
        """Returns iterable for leaf-nodes accessible from this node."""
        if not self.is_leaf:
            for child in self.children:
                yield from child.get_leaf_nodes()
        else:
            yield self

    def get_leaf_names(self) -> Generator[str | None, None, None]:
        for node in self.get_leaf_nodes():
            yield node.name

    def reroot_on_taxa(self, taxa: Iterable[str]) -> Newick:
        """Returns the Newick tree from this node, but rooted on the midpoint
        of the branch leading to one or more taxa. Note that the taxa are not
        required to form a clade. If the taxa do not form a monophyletic clade,
        then the outgroup will include more taxa than those passed to the
        function."""
        return _NewickGraph(self).reroot_on_taxa(taxa)

    def reroot_on_midpoint(self) -> Newick:
        """Returns the newick tree from this node, but rooted on the midpoint
        of the tree. That is to say that a root node is added at the exact
        midpoint of the longest path in the unrooted tree. If this midpoint
        lies at an existing internal node, then this node is made the root.

        Note that the sorting of nodes is not preserved, and that any
        uninformative nodes (lacking name/length, while connecting two
        other nodes, e.g. the old root) are spliced out.

        All nodes must have a length of zero or greater (no missing values
        are allowed), but note that rerooting behavior around nodes with
        length zero may yield unexpected results."""
        if len(list(self.get_leaf_nodes())) < 2:
            return self  # No meaningful way to reroot such trees

        return _NewickGraph(self).reroot_on_midpoint()

    def add_support(
        self,
        bootstraps: Iterable[Newick],
        fmt: str = "{Support}",
    ) -> Newick:
        """Adds support values to the current tree, based on a set of trees containing
        the same taxa. It is assumed that the support trees represent unrooted or
        arbitarily rooted trees, and no weight is given to the rooted topology of these
        trees.

        The main tree should itself be rooted, and the the toplogy and ordering of this
        tree is preserved, with node-names updated using the formatting string 'fmt'.

        Formatting is carried out using str.format, with these fields:
          {Support}    -- The total number of trees in which a clade is supported.
          {Percentage} -- The percentage of trees in which a clade is supported (float).
          {Fraction}   -- The fraction of trees in which a clade is supported (float).

        For example, typical percentage support-values can be realized by setting 'fmt'
        to the value "{Percentage:.0f}" to produce integer values.
        """
        clade_counts: dict[frozenset[str | None], int] = {}
        leaf_names_lst = list(self.get_leaf_names())
        leaf_names = frozenset(leaf_names_lst)
        if len(leaf_names) != len(leaf_names_lst):
            raise NewickError(
                "Cannot add support values to trees with duplicate leaf names"
            )

        n_bootstraps = 0
        for support_tree in bootstraps:
            n_bootstraps += 1
            support_tree_names = frozenset(support_tree.get_leaf_names())
            if leaf_names != support_tree_names:
                raise NewickError(
                    "Support tree does not contain same set of leaf nodes"
                )

            support_graph = _NewickGraph(support_tree)
            for clade in support_graph.get_clade_names():
                clade_counts[clade] = clade_counts.get(clade, 0) + 1

        return self._add_support(self, n_bootstraps, clade_counts, fmt)

    @classmethod
    def from_string(cls, string: str) -> Newick:
        """Parses a Newick string and returns a representation of the tree.
        See e.g. http://en.wikipedia.org/wiki/Newick_format

        Note that implicit nodes, such as (), (A,), and the like are not
        allowed, as they cannot always be represented/parsed in an unambigious
        manner. Thus all leaf nodes must have a name and/or a length."""
        tokens = _tokenize(string)
        if tokens and tokens[0] == "(":
            top_node = _parse_tokens(tokens)
        else:
            top_node = _parse_child(tokens)

        if tokens != [";"]:
            raise NewickParseError("Missing terminating semi-colon")

        return top_node

    def __lt__(self, other: object) -> bool:
        """See TotallyOrdered"""
        if not isinstance(other, Newick):
            return NotImplemented

        node_1 = (-self._weight, self.name or "", self.length or 0, self.children)
        node_2 = (-other._weight, other.name or "", other.length or 0, other.children)

        return node_1 < node_2

    def __hash__(self) -> int:
        """Hashing function, see 'hash'."""
        return self._hash

    def __repr__(self) -> str:
        """Representation corresponds to the Newick string for the (sub)tree,
        which can be parsed by 'from_string'."""
        return f"{self._to_str()};"

    def _to_str(self) -> str:
        fields: list[str] = []
        if self.children:
            fields.append("(")
            for child in self.children:
                fields.append(child._to_str())  # noqa: SLF001
                fields.append(",")
            fields.pop()
            fields.append(")")
        if self.name is not None:
            fields.append(str(self.name))
        if self.length is not None:
            fields.append(":")
            fields.append(str(self.length))
        return "".join(fields)

    def _add_support(
        self,
        node: Newick,
        total: int,
        clade_counts: dict[frozenset[str | None], int],
        fmt: str,
    ) -> Newick:
        """Recursively annotates a subtree with support values,
        excepting leaf nodes (where the name is preserved) and
        the root node (where the name is cleared)."""
        if node.is_leaf:
            return node

        clade = frozenset(leaf.name for leaf in node.get_leaf_nodes() if leaf.name)
        support = clade_counts.get(clade, 0)
        name = fmt.format(
            Support=support,
            Percentage=(support * 100.0) / (total or 1),
            Fraction=(support * 1.0) / (total or 1),
        )

        children: list[Newick] = []
        for child in node.children:
            children.append(self._add_support(child, total, clade_counts, fmt))

        return Newick(
            name=(None if (node is self) else name),
            length=node.length,
            children=children,
        )


################################################################################
################################################################################
# Functions related to NEWICK parsing

_TOKENIZER = re.compile("([():,;])")


def _tokenize(string: str) -> list[str]:
    result: list[str] = []
    for field in _TOKENIZER.split(string):
        field = field.strip()
        if field:
            result.append(field)
    return result


def _parse_tokens(tokens: list[str]) -> Newick:
    tokens.pop(0)
    child: Newick | None = None
    children: list[Newick] = []
    while tokens and tokens[0] not in ");":
        if tokens[0] == ",":
            if child is None:
                raise NewickParseError("Leaf nodes must have a name or a length")
            children.append(child)
            tokens.pop(0)
        child = _parse_child(tokens)

    if child is None:
        raise NewickParseError("Leaf nodes must have a name or a length")
    children.append(child)

    if not tokens or tokens[0] != ")":
        raise NewickParseError(
            "Malformed Newick string contains unbalanced parantheses"
        )
    tokens.pop(0)

    return _parse_child(tokens, children=children)


def _parse_child(tokens: list[str], children: Sequence[Newick] = ()) -> Newick:
    if tokens and tokens[0] == "(":
        return _parse_tokens(tokens)

    name: str | None = None
    length: str | None = None
    while tokens and (tokens[0] not in ",);"):
        if tokens[0] == ":":
            if length is not None:
                raise NewickParseError("Node has multiple length values")
            tokens.pop(0)
            if tokens[0] in ",);":
                raise NewickParseError("Missing length value")
            length = tokens.pop(0).strip()
        else:
            name = tokens.pop(0).strip()

    if not (name or length or children):
        raise NewickParseError("Parsing of implicit nodes not supported")

    return Newick(name=name, length=length, children=children)


################################################################################
################################################################################
# Class related to tree manipulations


class _NewickGraph:
    """Internal representation of an unrooted graph, allowing various forms of
    manipulation directly on the graph. To ensure that all manipulations can be
    carried out, it is required that branch-lengths are present for ALL branches,
    or for NO branches.

    Note that neither the root-length, nor node-ordering is preserved."""

    def __init__(self, node: Newick) -> None:
        self.names: dict[NodeID, NodeName] = {}
        self.connections: dict[NodeID, dict[NodeID, str | float | None]] = {}
        self.has_branch_lengths: bool | None = None
        self._collect_names_and_blengths(node)
        self.prune_uninformative_nodes()

    def _collect_names_and_blengths(self, c_node: Newick) -> None:
        c_node_id = id(c_node)

        self.set_name(c_node_id, c_node.name)
        for child in c_node.children:
            child_id = id(child)
            self.add_connection(c_node_id, child_id, child.length)
            self._collect_names_and_blengths(child)

    def is_leaf(self, node: NodeID) -> bool:
        """Returns true if the node is a leaf, defined as having a single connection."""
        return len(self.connections[node]) == 1

    def get_path_length(self, *nodes: NodeID) -> float | None:
        """Returns the length of a path through the graph. Calling the function
        with two nodes is the equivalent of getting the branch-length between
        those two nodes."""
        if not self.has_branch_lengths:
            return None

        path_length = 0.0
        for node_a, node_b in zip(nodes, nodes[1:]):
            segment_length = self.connections[node_a][node_b]
            if segment_length is None:
                raise AssertionError("segment_length is unexpectedly None")
            path_length += float(segment_length)

        return path_length

    def set_name(self, node_id: int, name: NodeName) -> None:
        self.names[node_id] = name

    def add_connection(
        self,
        node_id_a: NodeID,
        node_id_b: NodeID,
        blength: str | float | None = None,
    ) -> None:
        if (blength is not None) and float(blength) < 0:
            raise NewickError("Branch-lengths must be non-negative")
        elif (blength is not None) != self.has_branch_lengths:
            if self.has_branch_lengths is not None:
                raise NewickError("Tree contains branches with and without lengths")
            self.has_branch_lengths = blength is not None

        self.connections.setdefault(node_id_a, {})[node_id_b] = blength
        self.connections.setdefault(node_id_b, {})[node_id_a] = blength

    def remove_connection(
        self,
        node_a: NodeID,
        node_b: NodeID,
    ) -> str | float | None:
        length_a = self.connections[node_a].pop(node_b)
        length_b = self.connections[node_b].pop(node_a)
        assert length_a == length_b, (length_a, length_b)
        return length_a

    def remove_node(self, node: NodeID) -> None:
        connections = self.connections.pop(node)
        for node_b in connections:
            self.connections[node_b].pop(node)
        self.names.pop(node)

    def rebuild_tree(self, parent_id: NodeID, node_id: NodeID) -> Newick:
        """Rebuilds a newick tree starting at a node with id 'node_id' and a parent
        with id 'parent_id' (or the same value as 'node_id' if a root node).
        """
        children: list[Newick] = []
        for child_id in self.connections[node_id]:
            if child_id != parent_id:
                children.append(self.rebuild_tree(node_id, child_id))
        children.sort()

        blength = self.connections[parent_id].get(node_id)
        if isinstance(blength, float):
            blength = repr(blength)

        return Newick(name=self.names.get(node_id), length=blength, children=children)

    def prune_uninformative_nodes(self) -> None:
        """Removes nodes without names, and which are connected
        to two other nodes, extending the branch lengths of the
        two connected nodes. This process is repreated, until no
        further nodes are pruned. A rooted tree will typically
        contain just 1 such node, namely the old root node.

        For example, the tree "(A:5,(B:6):3);" would be reduced to
        the tree "(A:5,B:9);", whereas the trees "(A:5,(B:6)C:3);"
        and "(A:5,(B:6,C:2):3);" would not be pruned.

        For a node to be pruned, both adjacent nodes must have a
        length specified, or both must not have a length specified."""
        while True:
            for cur_node, connections in self.connections.items():
                if not self.names[cur_node] and (len(connections) == 2):
                    conn_a, conn_b = connections

                    blength = self.get_path_length(conn_a, cur_node, conn_b)

                    # Splice out the current node
                    self.remove_node(cur_node)
                    self.add_connection(conn_a, conn_b, blength)
                    break
            else:
                # Nothing was pruned this round, terminate
                break

    ################################################################################
    ################################################################################
    # Functions relating to NEWICK rooting on midpoint

    def reroot_on_midpoint(self) -> Newick:
        if not self.has_branch_lengths:
            raise NewickError(
                "Cannot reroot on midpoint for tree without branch-lengths"
            )

        longest_path, length = self._find_longest_path()
        root = self._create_root_at(longest_path, length / 2.0)

        return self.rebuild_tree(root, root)

    def _find_longest_path(self) -> tuple[list[NodeID], float]:
        """This function determines the longest non-overlapping path possible,
        and returns a list of the sequence of nodes in this path, as well as
        the total length of this path."""
        path_blengths: dict[frozenset[NodeID], float] = {}
        path_guides: dict[frozenset[NodeID], list[NodeID]] = {}

        def _collect_paths(
            guide: list[NodeID],
            length: float,
            p_node: NodeID,
            c_node: NodeID,
        ) -> None:
            path_len = self.get_path_length(p_node, c_node)
            assert path_len is not None

            length += path_len

            guide.append(c_node)
            key = frozenset(guide)
            path_blengths[key] = length
            path_guides[key] = guide

            for other in self.connections[c_node]:
                if other not in key:
                    _collect_paths(list(guide), length, c_node, other)

        for p_node, connections in self.connections.items():
            for c_node in connections:
                _collect_paths([p_node], 0, p_node, c_node)

        key, length = max(path_blengths.items(), key=lambda item: item[1])
        return path_guides[key], length

    def _create_root_at(
        self,
        path: list[NodeID],
        root_at: float,
    ) -> NodeID:
        """Finds the midpoint of a path through a tree, and
        either creates a new node at that point, or selects
        the node already present at that point (if any). The
        mid-point is assumed to be at distance of 'root_at'
        from the starting node.

        E.g. if the path is the longest path, and 'root_at' is
        half the length of this path, then this corresponds to
        rooting at the midpoint.

        The id of the new / selected node is returned. New
        nodes (if created) are always given the id None."""
        for c_node, n_node in zip(path, path[1:]):
            branch_length = self.get_path_length(c_node, n_node)
            assert branch_length is not None

            if branch_length > root_at:
                left_len = root_at
                right_len = branch_length - root_at

                self.remove_connection(c_node, n_node)
                self.add_connection(None, c_node, left_len)
                self.add_connection(None, n_node, right_len)

                return None
            elif branch_length == root_at:
                return n_node
            root_at -= branch_length

        raise AssertionError("midpoint not found")  # pragma: no coverage

    ################################################################################
    ################################################################################
    # Functions relating to NEWICK rooting on taxa

    def reroot_on_taxa(self, taxa: Iterable[str]) -> Newick:
        clades = self._collect_clades()
        root_on = self._collect_nodes_from_names(taxa)
        if not taxa:
            raise ValueError("No taxa in outgroup")

        # Because None is the id of the root atm:
        root = self._create_root_with_clade(clades, root_on)

        return self.rebuild_tree(root, root)

    def _collect_nodes_from_names(self, taxa: Iterable[str]) -> frozenset[NodeID]:
        nodes_by_names: dict[NodeName, list[NodeID]] = {}
        for node_id, name in self.names.items():
            if self.is_leaf(node_id):
                nodes_by_names.setdefault(name, []).append(node_id)

        selection: list[NodeID] = []
        for name in taxa:
            selection.extend(nodes_by_names.pop(name))

        if not nodes_by_names:
            raise ValueError("Cannot root on every taxa in tree")

        return frozenset(selection)

    def _collect_clades(self) -> dict[NodeID, dict[NodeID, frozenset[NodeID]]]:
        clades: dict[NodeID, dict[NodeID, frozenset[NodeID]]] = {}
        for node_a, connections in self.connections.items():
            for node_b in connections:
                self._collect_clade_from(clades, node_a, node_b)
        return clades

    def _collect_clade_from(
        self,
        cache: dict[NodeID, dict[NodeID, frozenset[NodeID]]],
        p_node: NodeID,
        c_node: NodeID,
    ) -> frozenset[NodeID]:
        c_clade = cache.setdefault(p_node, {}).get(c_node)
        if c_clade is not None:
            return c_clade

        clade: set[NodeID] = set()
        if self.is_leaf(c_node):
            clade.add(c_node)

        for n_node in self.connections[c_node]:
            if n_node != p_node:
                clade.update(self._collect_clade_from(cache, c_node, n_node))

        c_clade = frozenset(clade)
        cache[p_node][c_node] = c_clade
        return c_clade

    def _create_root_with_clade(
        self,
        clades: dict[NodeID, dict[NodeID, frozenset[NodeID]]],
        taxa: frozenset[NodeID],
    ) -> NodeID:
        root_key: tuple[NodeID, NodeID] | None = None
        root_clade: frozenset[NodeID] | None = None
        root_length: float | None = None
        for p_node, connections in clades.items():
            for n_node, clade in connections.items():
                if (
                    root_clade is None or len(clade) < len(root_clade)
                ) and taxa.issubset(clade):
                    root_key = (p_node, n_node)
                    root_clade = clade
                    root_length = self.get_path_length(p_node, n_node)

        assert root_key is not None
        p_node, n_node = root_key
        if root_length is not None:
            root_length = float(root_length) / 2.0

        self.remove_connection(p_node, n_node)
        self.add_connection(None, p_node, root_length)
        self.add_connection(None, n_node, root_length)

        return None

    ################################################################################
    ################################################################################
    # Functions relating to calculating bootstrap support

    def get_clade_names(self) -> set[frozenset[NodeName]]:
        result: set[frozenset[NodeName]] = set()
        for connections in self._collect_clades().values():
            for clade in connections.values():
                result.add(frozenset(self.names[node_id] for node_id in clade))
        return result
