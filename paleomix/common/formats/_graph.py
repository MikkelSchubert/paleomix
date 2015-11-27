#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
"""


"""

from paleomix.common.utilities import \
     safe_coerce_to_frozenset, \
     get_in, \
     set_in

from paleomix.common.formats import FormatError


class GraphError(FormatError):
    pass



class _Graph:
    """Internal representation of an unrooted graph, allowing various forms of
    manipulation directly on the graph. To ensure that all manipulations can be
    carried out, it is required that branch-lengths are present for ALL branches,
    or for NO branches.

    Note that neither the root-length, nor node-ordering is preserved."""

    def __init__(self):
        self.names        = {}
        self.connections  = {}
        self.has_branch_lengths = None


    def is_leaf(self, node):
        """Returns true if the node is a leaf, defined as having a single connection."""
        return len(self.connections[node]) == 1


    def get_path_length(self, *nodes):
        """Returns the length of a path through the graph. Calling the function
        with two nodes is the equivalent of getting the branch-length between
        those two nodes."""
        if not self.has_branch_lengths:
            return None

        path_length = 0.0
        for (node_a, node_b) in zip(nodes, nodes[1:]):
            segment_length = float(self.connections[node_a][node_b])
            path_length += segment_length

        return path_length


    def set_name(self, node_id, name):
        self.names[node_id] = name


    def add_connection(self, node_id_a, node_id_b, blength = None):
        if (blength is not None) and float(blength) < 0:
            raise GraphError("Branch-lengths must be non-negative")
        elif (blength is not None) != self.has_branch_lengths:
            if not self.has_branch_lengths is None:
                raise GraphError("Tree contains branches with and without lengths")
            self.has_branch_lengths = (blength is not None)

        set_in(self.connections, (node_id_a, node_id_b), blength)
        set_in(self.connections, (node_id_b, node_id_a), blength)


    def remove_connection(self, node_a, node_b):
        length_a = self.connections[node_a].pop(node_b)
        length_b = self.connections[node_b].pop(node_a)
        assert length_a == length_b, (length_a, length_b)
        return length_a


    def remove_node(self, node):
        connections = self.connections.pop(node)
        for node_b in connections:
            self.connections[node_b].pop(node)
        self.names.pop(node)


    def rebuild_tree(self, parent_id, node_id):
        """Rebuilds a tree starting at a node with id
        'node_id' and a parent with id 'parent_id' (or the
        same value as 'node_id' if a root node)."""
        raise NotImplementedError("Subclasses must implement 'rebuild_nodes'.") # pragma: no coverage


    def prune_uninformative_nodes(self):
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
            for (cur_node, connections) in self.connections.iteritems():
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
    ## Functions relating to NEWICK rooting on midpoint

    def reroot_on_midpoint(self):
        if not self.has_branch_lengths:
            raise GraphError("Cannot reroot on midpoint for tree without branch-lengths")

        longest_path, length = self._find_longest_path()
        root = self._create_root_at(longest_path, length / 2.0)

        return self.rebuild_tree(root, root)


    def _find_longest_path(self):
        """This function determines the longest non-overlapping path possible,
        and returns a list of the sequence of nodes in this path, as well as
        the total length of this path."""
        path_blengths = {}
        path_guides   = {}
        def _collect_paths(guide, length, p_node, c_node):
            length += self.get_path_length(p_node, c_node)

            guide.append(c_node)
            key                = frozenset(guide)
            path_blengths[key] = length
            path_guides[key ]  = guide

            for other in self.connections[c_node]:
                if other not in key:
                    _collect_paths(list(guide), length, c_node, other)

        for (p_node, connections) in self.connections.iteritems():
            for c_node in connections:
                _collect_paths([p_node], 0, p_node, c_node)

        key, length = max(path_blengths.iteritems(), key = lambda item: item[1])
        return path_guides[key], length


    def _create_root_at(self, path, root_at):
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
            branch_length = self.get_path_length(c_node, n_node)

            if (branch_length > root_at):
                left_len  = root_at
                right_len = branch_length - root_at

                self.remove_connection(c_node, n_node)
                self.add_connection(None, c_node, left_len)
                self.add_connection(None, n_node, right_len)

                return None
            elif branch_length == root_at:
                return n_node
            root_at -= branch_length
        assert False # pragma: no coverage


    ################################################################################
    ################################################################################
    ## Functions relating to NEWICK rooting on taxa

    def reroot_on_taxa(self, taxa):
        taxa  = safe_coerce_to_frozenset(taxa)
        if not taxa:
            raise ValueError("No taxa in outgroup")

        clades   = self._collect_clades()
        root_on  = self._collect_nodes_from_names(taxa)
        # Because None is the id of the root atm: # pylint: disable=W1111
        root     = self._create_root_with_clade(clades, root_on)

        return self.rebuild_tree(root, root)


    def _collect_nodes_from_names(self, taxa):
        known_taxa = set()
        for (node_id, name) in self.names.iteritems():
            if self.is_leaf(node_id):
                known_taxa.add(name)

        unknown_taxa = taxa - known_taxa
        if unknown_taxa:
            raise ValueError("Cannot root on unknown taxa: %s" % (", ".join(unknown_taxa),))
        elif not (known_taxa - taxa):
            raise ValueError("Cannot root on every taxa in tree")

        return frozenset(key for (key, name) in self.names.iteritems() if name in taxa)


    def _collect_clades(self):
        clades = {}
        for (node_a, connections) in self.connections.iteritems():
            for node_b in connections:
                self._collect_clade_from(clades, node_a, node_b)
        return clades


    def _collect_clade_from(self, cache, p_node, c_node):
        c_clade = get_in(cache, (p_node, c_node), set())
        if not c_clade:
            if self.is_leaf(c_node):
                c_clade.add(c_node)

            for n_node in self.connections[c_node]:
                if n_node != p_node:
                    c_clade.update(self._collect_clade_from(cache, c_node, n_node))
            set_in(cache, (p_node, c_node), frozenset(c_clade))
        return c_clade


    def _create_root_with_clade(self, clades, taxa):
        root_key, root_clade, root_length = None, None, None
        for (p_node, connections) in clades.iteritems():
            for (n_node, clade) in connections.iteritems():
                if (root_clade is None) or (len(clade) < len(root_clade)):
                    if taxa.issubset(clade):
                        root_key    = (p_node, n_node)
                        root_clade  = clade
                        root_length = self.get_path_length(p_node, n_node)

        p_node, n_node = root_key
        if root_length is not None:
            root_length = float(root_length) / 2.0

        self.remove_connection(p_node, n_node)
        self.add_connection(None, p_node, root_length)
        self.add_connection(None, n_node, root_length)

        return None


    ################################################################################
    ################################################################################
    ## Functions relating to calculating bootstrap support
    def get_clade_names(self):
        result = set()
        for (_, connections) in self._collect_clades().iteritems():
            for (_, clade) in connections.iteritems():
                result.add(frozenset(self.names[node_id] for node_id in clade))
        return result
