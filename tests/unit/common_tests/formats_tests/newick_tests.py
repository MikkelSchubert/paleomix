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
from nose.tools import \
     assert_equal, \
     assert_raises

from pypeline.common.formats.newick import \
     Newick, \
     NewickError, \
     NewickParseError

from pypeline.common.testing import \
     assert_list_equal


############################################################################
############################################################################
## Constructor

def test_newick__constructor__name():
    node = Newick(name = "AbC")
    assert_equal(node.name, "AbC")

def test_newick__constructor__children_set_in_internal_nodes():
    node = Newick(name = "Leaf")
    top_node = Newick(children = [node])
    assert_equal(top_node.children, [node])

def test_newick__constructor__children_not_set_in_leafs():
    node = Newick(name = "Leaf")
    assert_equal(node.children, None)

def test_newick__constructor__is_leaf_true_for_leafs():
    node = Newick(name = "Another Leaf")
    assert node.is_leaf

def test_newick__constructor__is_leaf_false_for_internal_nodes():
    node = Newick(name = "Leaf")
    top_node = Newick(children = [node])
    assert not top_node.is_leaf

def test_newick__constuctor__leafs_must_have_name_or_length():
    assert_raises(NewickError, Newick, children = None)

def test_newick__constructor__internal_nodes_must_have_children():
    assert_raises(NewickError, Newick, children = [])


############################################################################
############################################################################
## get_leafs

def test_newick__from_leafs__leaf_returns_self():
    node = Newick(name = "Leaf")
    assert_list_equal(node.get_leafs(), [node])

def test_newick__from_leafs__internal_node_returns_leafs():
    node_a = Newick(name = "Leaf A")
    node_b = Newick(name = "Leaf B")
    top_node = Newick(children = [node_a, node_b])
    assert_list_equal(top_node.get_leafs(), [node_a, node_b])

def test_newick__from_leafs__complex_case():
    node_a = Newick(name = "Leaf A")
    node_b = Newick(name = "Leaf B")
    node_c = Newick(name = "Leaf C")
    sub_a  = Newick(children = [node_b, node_c])
    top_node = Newick(children = [node_a, sub_a])
    assert_list_equal(top_node.get_leafs(), [node_a, node_b, node_c])


############################################################################
############################################################################
## from_string

def test_newick__parse__minimal_newick__name_only():
    top_node = Newick(name = "A")
    assert_equal(Newick.from_string("A;"), top_node)

def test_newick__parse__single_taxa():
    child_node = Newick(name = "Ab")
    top_node   = Newick(children = [child_node])
    assert_equal(Newick.from_string("(Ab);"), top_node)

def test_newick__parse__two_taxa():
    child_node_1 = Newick(name = "A")
    child_node_2 = Newick(name = "Bc")
    top_node   = Newick(children = [child_node_1, child_node_2])
    assert_equal(Newick.from_string("(A,Bc);"), top_node)

def test_newick__parse__three_taxa():
    child_node_1 = Newick(name = "A")
    child_node_2 = Newick(name = "Bc")
    child_node_3 = Newick(name = "DeF")
    top_node   = Newick(children = [child_node_1, child_node_2, child_node_3])
    assert_equal(Newick.from_string("(A,Bc,DeF);"), top_node)

def test_newick__parse__ignore_whitespace():
    assert_equal(Newick.from_string("(A,B);"), Newick.from_string("(A, B);"))

def test_newick__parse__missing_semicolon():
    assert_raises(NewickParseError, Newick.from_string, "()")


def test_newick__parse__subnode__single_taxa():
    child_node_1  = Newick(name = "A")
    child_node_2a = Newick(name = "B")
    child_node_2  = Newick(children = [child_node_2a])
    top_node   = Newick(children = [child_node_1, child_node_2])
    assert_equal(Newick.from_string("(A,(B));"), top_node)

def test_newick__parse__subnode__two_taxa():
    child_node_1  = Newick(name = "A")
    child_node_2a = Newick(name = "B")
    child_node_2b = Newick(name = "C")
    child_node_2  = Newick(children = [child_node_2a, child_node_2b])
    top_node   = Newick(children = [child_node_1, child_node_2])
    assert_equal(Newick.from_string("(A,(B,C));"), top_node)



############################################################################
############################################################################
## Malformed newick strings

def test_newick__malformed__unbalanced_parantheses():
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C);")

def test_newick__malformed__mismatched_parantheses():
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C();")

def test_newick__malformed__missing_parantheses():
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C))")

def test_newick__malformed__missing_length():
    assert_raises(NewickParseError, Newick.from_string, "(A:,(B,C));")
    assert_raises(NewickParseError, Newick.from_string, "(A,(B:,C));")
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C:));")
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C):);")
    assert_raises(NewickParseError, Newick.from_string, "(A,(B,C)):;")

def test_newick__malformed__multiple_lengths():
    assert_raises(NewickParseError, Newick.from_string, "(A:1:2);")


############################################################################
############################################################################
## Implicit leafs are not supported (due to problems with ambiguiety)

def test_newick__parse__first_taxa_unnamed():
    assert_raises(NewickError, Newick.from_string, "(,A);")

def test_newick__parse__second_taxa_unnamed():
    assert_raises(NewickError, Newick.from_string, "(A,);")

def test_newick__parse__two_taxa_unnamed():
    assert_raises(NewickError, Newick.from_string, "(,);")

def test_newick__parse__three_taxa_unnamed():
    assert_raises(NewickError, Newick.from_string, "(,,);")


############################################################################
############################################################################
## Empty non-leaf nodes are not allowed, as their interpretation is unclear

def test_newick__parse__minimal_newick__implicit_nodes():
    assert_raises(NewickParseError, Newick.from_string, "();")

def test_newick__parse__subnode__empty():
    assert_raises(NewickParseError, Newick.from_string, "(A,());")


############################################################################
############################################################################
# The following tests are derived from the wikipedia description of the
# newick format: http://en.wikipedia.org/wiki/Newick_format#Examples

def test_newick__wikipedia_example_1():
    # no nodes are named, this format is not supported here!
    assert_raises(NewickError, Newick.from_string, "(,,(,));")

def test_newick__wikipedia_example_2():
    # leaf nodes are named
    taxa_d = Newick(name = "D")
    taxa_c = Newick(name = "C")
    taxa_sub = Newick(children = [taxa_c, taxa_d])
    taxa_b = Newick(name = "B")
    taxa_a = Newick(name = "A")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub])
    assert_equal(Newick.from_string("(A,B,(C,D));"), top_node)

def test_newick__wikipedia_example_3():
    # all nodes are named
    taxa_d = Newick(name = "D")
    taxa_c = Newick(name = "C")
    taxa_sub = Newick(children = [taxa_c, taxa_d], name = "E")
    taxa_b = Newick(name = "B")
    taxa_a = Newick(name = "A")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub], name = "F")
    assert_equal(Newick.from_string("(A,B,(C,D)E)F;"), top_node)

def test_newick__wikipedia_example_4():
    # all but root node have a distance to parent
    taxa_d = Newick(length = "0.4")
    taxa_c = Newick(length = "0.3")
    taxa_sub = Newick(children = [taxa_c, taxa_d], length = "0.5")
    taxa_b = Newick(length = "0.2")
    taxa_a = Newick(length = "0.1")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub])
    assert_equal(Newick.from_string("(:0.1,:0.2,(:0.3,:0.4):0.5);"), top_node)

def test_newick__wikipedia_example_5():
    # all have a distance to parent
    taxa_d = Newick(length = "0.4")
    taxa_c = Newick(length = "0.3")
    taxa_sub = Newick(children = [taxa_c, taxa_d], length = "0.5")
    taxa_b = Newick(length = "0.2")
    taxa_a = Newick(length = "0.1")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub], length = "0.0")
    assert_equal(Newick.from_string("(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"), top_node)

def test_newick__wikipedia_example_6():
    # distances and leaf names (popular)
    taxa_d = Newick(length = "0.4", name = "D")
    taxa_c = Newick(length = "0.3", name = "C")
    taxa_sub = Newick(children = [taxa_c, taxa_d], length = "0.5")
    taxa_b = Newick(length = "0.2", name = "B")
    taxa_a = Newick(length = "0.1", name = "A")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub])
    assert_equal(Newick.from_string("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"), top_node)

def test_newick__wikipedia_example_7():
    # distances and all names
    taxa_d = Newick(length = "0.4", name = "D")
    taxa_c = Newick(length = "0.3", name = "C")
    taxa_sub = Newick(children = [taxa_c, taxa_d], length = "0.5", name = "E")
    taxa_b = Newick(length = "0.2", name = "B")
    taxa_a = Newick(length = "0.1", name = "A")
    top_node = Newick(children = [taxa_a, taxa_b, taxa_sub], name = "F")
    assert_equal(Newick.from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"), top_node)

def test_newick__wikipedia_example_8():
    # a tree rooted on a leaf node (rare)
    taxa_b = Newick(length = "0.2", name = "B")
    taxa_c = Newick(length = "0.3", name = "C")
    taxa_d = Newick(length = "0.4", name = "D")
    node_e = Newick(length = "0.5", name = "E", children = [taxa_c, taxa_d])
    node_f = Newick(length = "0.1", name = "F", children = [taxa_b, node_e])
    node_a = Newick(name = "A", children = [node_f])
    assert_equal(Newick.from_string("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"), node_a)



############################################################################
############################################################################
## str / repr

def test_newick__str__repr_equal_to_str():
    node_a = Newick(name = "A", length = "123")
    node_b = Newick(name = "B", length = "4567")
    top_node = Newick(children = [node_a, node_b])
    assert_equal(str(top_node), "(A:123,B:4567);")

def test_newick__str__single_leaf_should_not_be_followed_by_comma():
    node = Newick(name = "A")
    top_node = Newick(children = [node])
    assert_equal(str(top_node), "(A);")

def test_newick__wikipedia_examples__str_equality():
    def _newick_str_input_equals_output(nwk_str):
        nodes  = Newick.from_string(nwk_str)
        result = str(nodes)
        assert_equal(result, nwk_str)

    # 2. leaf nodes are named
    yield _newick_str_input_equals_output, "(A,B,(C,D));"
    # 3. all nodes are named
    yield _newick_str_input_equals_output, "(A,B,(C,D)E)F;"
    # 4. all but root node have a distance to parent
    yield _newick_str_input_equals_output, "(:0.1,:0.2,(:0.3,:0.4):0.5);"
    # 5. all have a distance to parent
    yield _newick_str_input_equals_output, "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"
    # 6. distances and leaf names (popular)
    yield _newick_str_input_equals_output, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    # 7. distances and all names
    yield _newick_str_input_equals_output, "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
