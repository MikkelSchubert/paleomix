#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
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
# Required due to use of NotImplementedError in setattr:
# pylint: disable=R0921
import re

from pypeline.common.utilities import set_in


class NewickError(RuntimeError):
    pass


class NewickParseError(NewickError):
    pass


class Newick(object):
    """Immutable object representing a Newick node.

    Nodes are classified as either internal nodes (have children),
    or leaf nodes (does not have children). A node MUST either have
    1 or more child-nodes, or have a name and/or a length. This is to
    ensure that nodes can be represented in an unambigious manner
    using the Newick format.

    No assumptions are made about the type of the 'name' and the 'length'
    properties, and these are simply converted into strings when the
    Newick string is generated.
    """

    def __init__(self, name = None, length = None, children = None):
        object.__init__(self)
        object.__setattr__(self, "name",     name)
        object.__setattr__(self, "length",   length)
        object.__setattr__(self, "children", tuple(children or ()))

        # Ensure that these values are hashable
        hash(self.name)
        hash(self.length)

        for child in self.children:
            if not isinstance(child, Newick):
                raise TypeError("Child nodes must be Newick nodes")

        if not (self.children or self.name or self.length):
            raise NewickError("Leaf nodes MUST have either a name or a length")


    @property
    def is_leaf(self):
        """Returns true if the node is a leaf (has no children)."""
        return not self.children


    def get_leaf_nodes(self):
        """Returns iterable for leaf-nodes accessible from this node."""
        if not self.is_leaf:
            for child in self.children:
                for leaf in child.get_leaf_nodes():
                    yield leaf
        else:
            yield self


    def reroot_on_midpoint(self):
        """Returns the newick tree from this node, but rooted on the midpoint
        of the tree. That is to say that a root node is added at the exact
        midpoint of the longest path in the unrooted tree. If this midpoint
        lies at an existing internal node, then this node is made the root.

        Note that the sorting of nodes is not preserved, and that any
        uninformative nodes (lacking name/length, while connecting two
        other nodes, e.g. the old root) are spliced out."""
        return _reroot_on_midpoint(self)


    @classmethod
    def from_string(cls, string):
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


    def __cmp__(self, other):
        return cmp((self.name,
                    self.length,
                    self.children),
                   (other.name,
                    other.length,
                    other.children))


    def __hash__(self):
        return hash((self.name,
                    self.length,
                    self.children))


    def __repr__(self):
        return "%s;" % (self._to_str(),)


    def _to_str(self):
        fields = []
        if self.children:
            fields.append("(")
            for child in self.children:
                fields.append(child._to_str()) # pylint: disable=W0212
                fields.append(",")
            fields.pop()
            fields.append(")")
        if self.name is not None:
            fields.append(str(self.name))
        if self.length is not None:
            fields.append(":")
            fields.append(str(self.length))
        return "".join(fields)


    def __setattr__(self, _name, _value):
        raise NotImplementedError("Newick nodes are immutable")


    def __delattr__(self, _name):
        raise NotImplementedError("Newick nodes are immutable")



################################################################################
################################################################################
## Functions relating to NEWICK rooting

def _reroot_on_midpoint(node):
    """See Newick.reroot_at_midpoint."""
    # Dictionary of object IDs of nodes to node names
    names   = _collect_node_names(node, {})
    # Dictionary of object IDs to object IDs to branch lengths
    blengths = _collect_branch_lengths(node,{})
    blengths = _prune_uninformative_nodes(blengths, names)
    longest_path, length = _find_longest_path(blengths)
    root = _create_root_at(blengths, longest_path, length / 2.0)

    return _rebuild_nodes(blengths, names, root, root)


def _collect_node_names(node, table):
    """Returns a dictionary of id(node) -> node.name."""
    table[id(node)] = node.name
    for child in node.children:
        _collect_node_names(child, table)
    return table


def _collect_branch_lengths(node, table):
    """Returns a dictionary
      { id(node_a) : { id(node_b) : branch length } }
    containing the branch lengths for all pairs of nodes where
      node_a is not node_b
    and one node is a child node of the other node."""
    node_id = id(node)
    for child in node.children:
        if child.length is None:
            raise ValueError("Branch-lengths must be specified for ALL nodes")

        child_id = id(child)
        length = float(child.length)

        set_in(table, (node_id, child_id), length)
        set_in(table, (child_id, node_id), length)

        _collect_branch_lengths(child, table)

    return table


def _prune_uninformative_nodes(blengths, names):
    """Removes nodes without names, and which are connected
    to two other nodes, extending the branch lengths of the
    two connected nodes. This process is repreated, until no
    further nodes are pruned. A rooted tree will typically
    contain just 1 such node, namely the old root node.

    For example, the tree "(A:5,(B:6):3);" would be reduced to
    the tree "(A:5,B:9);", whereas the trees "(A:5,(B:6)C:3);"
    and "(A:5,(B:6,C:2):3);" would not be pruned."""
    while True:
        for (cur_node, connections) in blengths.iteritems():
            if not names[cur_node] and (len(connections) == 2):
                other_node_a, other_node_b = connections

                length \
                  = blengths[other_node_a].pop(cur_node) \
                  + blengths[other_node_b].pop(cur_node)

                set_in(blengths, (other_node_a, other_node_b), length)
                set_in(blengths, (other_node_b, other_node_a), length)

                del blengths[cur_node]
                break
        else:
            # Nothing was pruned this round, terminate
            break

    return blengths


def _find_longest_path(blengths):
    """Given a dictionary of branch-lengths (see _collect_branch_lengths)
    this function determines the longest non-overlapping path possible,
    and returns a list of the sequence of nodes in this path, as well as
    the total length of this path."""
    path_blengths = {}
    path_guides   = {}
    def _collect_paths(guide, length, p_node, c_node):
        guide.append(c_node)
        key = frozenset(guide)

        length += blengths[p_node][c_node]
        path_blengths[key] = length
        path_guides[key]  = guide

        for other in blengths[c_node]:
            if other not in key:
                _collect_paths(list(guide), length, c_node, other)

    for (p_node, branches) in blengths.iteritems():
        for c_node in branches:
            _collect_paths([p_node], 0, p_node, c_node)

    key, length = max(path_blengths.iteritems(), key = lambda item: item[1])
    return path_guides[key], length


def _create_root_at(lengths, path, root_at):
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
    for (c_node, n_node) in zip(path, path[1:]):
        branch_length = lengths[c_node][n_node]

        if branch_length > root_at:
            del lengths[c_node][n_node]
            del lengths[n_node][c_node]

            left_len  = root_at
            right_len = branch_length - root_at

            set_in(lengths, (None, c_node), left_len)
            set_in(lengths, (c_node, None), left_len)
            set_in(lengths, (None, n_node), right_len)
            set_in(lengths, (n_node, None), right_len)

            return None
        elif branch_length == root_at:
            return n_node
        root_at -= branch_length

    assert False # pragma: no coverage


def _rebuild_nodes(blengths, names, parent_id, node_id):
    """Rebuilds a newick tree starting at a node with id
    'node_id' and a parent with id 'parent_id' (or the
    same value as 'node_id' if a root node).

    Takes a dict of branch-lengths as returned by the
    function _collect_branch_lengths, and a dict of node
    names as returned by _collect_node_names."""

    children = []
    for child_id in blengths[node_id]:
        if child_id != parent_id:
            children.append(_rebuild_nodes(blengths, names, node_id, child_id))
    children.sort()

    blength = blengths.get(parent_id).get(node_id)
    if blength is not None:
        blength = repr(blength)

    return Newick(name     = names.get(node_id),
                  length   = blength,
                  children = children)




################################################################################
################################################################################
## Functions related to NEWICK parsing

_TOKENIZER = re.compile("([():,;])")
_NODE_KEYS = frozenset(("name", "length", "children"))


def _tokenize(string):
    result = []
    for field in _TOKENIZER.split(string):
        field = field.strip()
        if field:
            result.append(field)
    return result


def _parse_tokens(tokens):
    assert tokens and tokens[0] == "("

    tokens.pop(0)
    child, children = None, []
    while tokens and (tokens[0] not in ");"):
        if tokens[0] == ",":
            children.append(child)
            tokens.pop(0)
        child = _parse_child(tokens)
    children.append(child)

    if any(child is None for child in children):
        raise NewickParseError("Implicit leaf nodes (no name OR length) are not allowed")
    elif not tokens or (tokens[0] != ")"):
        raise NewickParseError("Malformed Newick string, contains unbalanced parantheses")
    tokens.pop(0)

    return _parse_child(tokens, children = children)


def _parse_child(tokens, children = None):
    if tokens and tokens[0] == "(":
        return _parse_tokens(tokens)

    name, length = None, None
    while tokens and (tokens[0] not in ",);"):
        if (tokens[0] == ":"):
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

    return Newick(name     = name,
                  length   = length,
                  children = children)
