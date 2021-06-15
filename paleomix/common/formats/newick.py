#!/usr/bin/python3
#
# Copyright (c) 2012 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import re
from typing import (
    Any,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Optional,
    Tuple,
    Union,
)

from paleomix.common.formats._graph import Graph, GraphError
from paleomix.common.utilities import Immutable, TotallyOrdered

NodeID = Optional[int]
NodeName = Optional[str]


class NewickError(GraphError):
    pass


class NewickParseError(NewickError):
    """Exception raised if errors occur during parsing
    of Newick strings."""

    pass


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
    contraints apply when unrooting/rerooting trees (see below)."""

    __slots__ = ["name", "length", "children", "_hash", "_weight"]
    name: Optional[str]
    length: Optional[Union[float, str]]
    children: Tuple["Newick", ...]
    _hash: int
    _weight: int

    def __init__(
        self,
        name: Optional[str] = None,
        length: Optional[Union[str, float]] = None,
        children: Iterable["Newick"] = (),
    ):
        children = tuple(children)
        name = name or None
        length = length or None

        Immutable.__init__(
            self,
            name=name,
            length=length,
            children=children,
            _hash=hash((name, length, children)),
        )

        if not (self.children or self.name or self.length):
            raise NewickError("Leaf nodes MUST have either a name or a length")

        weight = 0
        for child in self.children:
            if not isinstance(child, Newick):
                raise TypeError("Child nodes must be Newick nodes")
            weight += 1
        object.__setattr__(self, "_weight", weight)

    @property
    def is_leaf(self):
        """Returns true if the node is a leaf (has no children)."""
        return not self.children

    def get_leaf_nodes(self) -> Iterator["Newick"]:
        """Returns iterable for leaf-nodes accessible from this node."""
        if not self.is_leaf:
            for child in self.children:
                for leaf in child.get_leaf_nodes():
                    yield leaf
        else:
            yield self

    def get_leaf_names(self):
        for node in self.get_leaf_nodes():
            yield node.name

    def reroot_on_taxa(self, taxa: Iterable[str]) -> "Newick":
        """Returns the Newick tree from this node, but rooted on the midpoint
        of the branch leading to one or more taxa. Note that the taxa are not
        required to form a clade. If the taxa do not form a monophyletic clade,
        then the outgroup will include more taxa than those passed to the
        function."""
        return _NewickGraph(self).reroot_on_taxa(taxa)

    def reroot_on_midpoint(self):
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

    def add_support(self, bootstraps: Iterable["Newick"], fmt: str = "{Support}"):
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
        clade_counts = {}  # type: Dict[FrozenSet[Optional[str]], int]
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
    def from_string(cls, string: str) -> "Newick":
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

    def __lt__(self, other: Any):
        """See TotallyOrdered"""
        if not isinstance(other, Newick):
            return NotImplemented

        return (-self._weight, self.name, self.length, self.children) < (
            -other._weight,
            other.name,
            other.length,
            other.children,
        )

    def __hash__(self):
        """Hashing function, see 'hash'."""
        return self._hash

    def __repr__(self):
        """Representation corresponds to the Newick string for the (sub)tree,
        which can be parsed by 'from_string'."""
        return "%s;" % (self._to_str(),)

    def _to_str(self):
        fields = []  # type: List[str]
        if self.children:
            fields.append("(")
            for child in self.children:
                fields.append(child._to_str())
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
        node: "Newick",
        total: int,
        clade_counts: Dict[FrozenSet[Optional[str]], int],
        fmt: str,
    ) -> "Newick":
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

        children = []  # type: List[Newick]
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


def _tokenize(string: str) -> List[str]:
    result = []  # type: List[str]
    for field in _TOKENIZER.split(string):
        field = field.strip()
        if field:
            result.append(field)
    return result


def _parse_tokens(tokens: List[str]) -> Newick:
    tokens.pop(0)
    child = None  # type: Optional[Newick]
    children = []  # type: List[Newick]
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


def _parse_child(tokens: List[str], children: List[Newick] = []):
    if tokens and tokens[0] == "(":
        return _parse_tokens(tokens)

    name: Optional[str] = None
    length: Optional[str] = None
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


class _NewickGraph(Graph):
    def __init__(self, node: Newick):
        Graph.__init__(self)
        self._collect_names_and_blengths(node)
        self.prune_uninformative_nodes()

    def _collect_names_and_blengths(self, c_node: Newick):
        c_node_id = id(c_node)

        self.set_name(c_node_id, c_node.name)
        for child in c_node.children:
            child_id = id(child)
            self.add_connection(c_node_id, child_id, child.length)
            self._collect_names_and_blengths(child)

    def rebuild_tree(self, parent_id: NodeID, node_id: NodeID) -> "Newick":
        """Rebuilds a newick tree starting at a node with id
        'node_id' and a parent with id 'parent_id' (or the
        same value as 'node_id' if a root node)."""

        children = []  # type: List[Newick]
        for child_id in self.connections[node_id]:
            if child_id != parent_id:
                children.append(self.rebuild_tree(node_id, child_id))
        children.sort()

        blength = self.connections[parent_id].get(node_id)
        if isinstance(blength, float):
            blength = repr(blength)

        return Newick(name=self.names.get(node_id), length=blength, children=children)
