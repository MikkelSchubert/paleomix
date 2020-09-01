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
import pytest

from paleomix.common.formats.newick import (
    GraphError,
    Newick,
    NewickError,
    NewickParseError,
)

###############################################################################
###############################################################################
# Constructor


def test_newick__constructor__name():
    node = Newick(name="AbC")
    assert node.name == "AbC"


def test_newick__constructor__children_set_in_internal_nodes():
    node = Newick(name="Leaf")
    top_node = Newick(children=[node])
    assert top_node.children == (node,)


def test_newick__constructor__children_not_set_in_leaf_nodes():
    node = Newick(name="Leaf")
    assert node.children == ()


def test_newick__constructor__is_leaf_true_for_leaf_nodes():
    node = Newick(name="Another Leaf")
    assert node.is_leaf


def test_newick__constructor__is_leaf_false_for_internal_nodes():
    node = Newick(name="Leaf")
    top_node = Newick(children=[node])
    assert not top_node.is_leaf


def test_newick__constuctor__leaf_nodes_must_have_name_or_length():
    with pytest.raises(NewickError):
        Newick(children=None)


def test_newick__constructor__internal_nodes_must_have_children():
    with pytest.raises(NewickError):
        Newick(children=[])


def test_newick__constructor__children_must_be_newick():
    with pytest.raises(TypeError):
        Newick(children=["A", "B"])


###############################################################################
###############################################################################
# get_leaf_nodes


def test_newick__get_leaf_nodes__leaf_returns_self():
    node = Newick(name="Leaf")
    assert list(node.get_leaf_nodes()) == list([node])


def test_newick__get_leaf_nodes__internal_node_returns_leaf_nodes():
    node_a = Newick(name="Leaf A")
    node_b = Newick(name="Leaf B")
    top_node = Newick(children=[node_a, node_b])
    assert list(top_node.get_leaf_nodes()) == list([node_a, node_b])


def test_newick__get_leaf_nodes__complex_case():
    node_a = Newick(name="Leaf A")
    node_b = Newick(name="Leaf B")
    node_c = Newick(name="Leaf C")
    sub_a = Newick(children=[node_b, node_c])
    top_node = Newick(children=[node_a, sub_a])
    assert list(top_node.get_leaf_nodes()) == list([node_a, node_b, node_c])


###############################################################################
###############################################################################
# get_leaf_nodes


def test_newick__get_leaf_names__leaf_returns_self():
    node = Newick(name="Leaf")
    assert list(node.get_leaf_names()) == list(["Leaf"])


def test_newick__get_leaf_names__internal_node_returns_leaf_nodes():
    node_a = Newick(name="Leaf A")
    node_b = Newick(name="Leaf B")
    top_node = Newick(children=[node_a, node_b])
    assert list(top_node.get_leaf_names()) == list(["Leaf A", "Leaf B"])


def test_newick__get_leaf_names__complex_case():
    node_a = Newick(name="Leaf A")
    node_b = Newick(name="Leaf B")
    node_c = Newick(name="Leaf C")
    sub_a = Newick(children=[node_b, node_c])
    top_node = Newick(children=[node_a, sub_a])
    assert list(top_node.get_leaf_names()) == list(["Leaf A", "Leaf B", "Leaf C"])


###############################################################################
###############################################################################
# reroot_on_taxa


def test_newick__reroot_on_taxa__single_taxa():
    source = Newick.from_string("((A,B),C);")
    expected = Newick.from_string("((B,C),A);")
    assert expected == source.reroot_on_taxa("A")


def test_newick__reroot_on_taxa__single_taxa_with_branch_lengths():
    source = Newick.from_string("((A:4,B:3):2,C:1);")
    expected = Newick.from_string("((B:3,C:3.0):2.0,A:2.0);")
    assert expected == source.reroot_on_taxa("A")


def test_newick__reroot_on_taxa__multiple_taxa__clade():
    source = Newick.from_string("((A,(B,C)),(D,E));")
    expected = Newick.from_string("(((D,E),A),(B,C));")
    assert expected == source.reroot_on_taxa(("B", "C"))


def test_newick__reroot_on_taxa__multiple_taxa__paraphylogeny():
    source = Newick.from_string("((B,C),((D,E),A));")
    expected = Newick.from_string("(((B,C),A),(D,E));")
    assert expected == source.reroot_on_taxa(("A", "C"))


def test_newick__reroot_on_taxa__no_taxa():
    source = Newick.from_string("((B,C),((D,E),A));")
    with pytest.raises(ValueError):
        source.reroot_on_taxa(())


def test_newick__reroot_on_taxa__unknown_taxa():
    source = Newick.from_string("((B,C),((D,E),A));")
    with pytest.raises(ValueError):
        source.reroot_on_taxa(("A", "Z"))


def test_newick__reroot_on_taxa__no_non_root_taxa():
    source = Newick.from_string("((B,C),((D,E),A));")
    with pytest.raises(ValueError):
        source.reroot_on_taxa(("A", "B", "C", "D", "E"))


###############################################################################
###############################################################################
# reroot_on_midpoint


def test_newick__reroot_on_midpoint__single_node():
    source = Newick.from_string("(A:3.0);")
    expected = Newick.from_string("(A:3.0);")
    assert expected == source.reroot_on_midpoint()


def test_newick__reroot_on_midpoint__two_nodes():
    source = Newick.from_string("(A:3.0,B:8.0);")
    rerooted = source.reroot_on_midpoint()
    expected = Newick.from_string("(A:5.5,B:5.5);")
    assert expected == rerooted


def test_newick__reroot_on_midpoint__two_clades():
    source = Newick.from_string("((A:7,B:2):1,(C:1,D:0.5):2);")
    rerooted = source.reroot_on_midpoint()
    expected = Newick.from_string("(((C:1,D:0.5):3.0,B:2):1.5,A:5.5);")
    assert expected == rerooted


def test_newick__reroot_on_midpoint__nested_clades():
    source = Newick.from_string("((A:2,(B:2,C:3):4):1,(D:1,E:0.5):2);")
    rerooted = source.reroot_on_midpoint()
    expected = Newick.from_string("(((D:1,E:0.5):3.0,A:2):1.5,(B:2,C:3):2.5);")
    assert expected == rerooted


def test_newick__reroot_on_midpoint__reroot_on_internal_node():
    source = Newick.from_string("((A:5.0,B:1.0)C:2.0,D:3.0);")
    rerooted = source.reroot_on_midpoint()
    expected = Newick.from_string("(A:5.0,B:1.0,D:5.0)C;")
    assert expected == rerooted


_INVALID_BRANCH_LENGTHS = (
    "(A,B);",  # No branch lengths
    "(A:7,B);",  # Length missing for leaf node
    "(A:7,(B:3));",  # Length missing for inner node
    "(A:7,(B:3):-1);",  # Negative branch length
    "(A:7,B:-1);",  # Negative leaf length
)


@pytest.mark.parametrize("newick", _INVALID_BRANCH_LENGTHS)
def test_newick__reroot_on_midpoint__invalid_branch_lengths(newick):
    source = Newick.from_string(newick)
    with pytest.raises(GraphError):
        source.reroot_on_midpoint()


###############################################################################
###############################################################################
# add_support


def test_newick__add_support__no_trees():
    main_tree = Newick.from_string("(((A,B),C),D);")
    expected = Newick.from_string("(((A,B)0,C)0,D);")
    result = main_tree.add_support([])
    assert expected == result


def test_newick__add_support__single_identical_tree():
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [Newick.from_string("(((A,B),C),D);")]
    expected = Newick.from_string("(((A,B)1,C)1,D);")
    result = main_tree.add_support(bootstraps)
    assert expected == result


def test_newick__add_support__single_identical_tree__different_rooting():
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [Newick.from_string("(((C,D),B),A);")]
    expected = Newick.from_string("(((A,B)1,C)1,D);")
    result = main_tree.add_support(bootstraps)
    assert expected == result


def test_newick__add_support__multiple_trees__different_topologies():
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [
        Newick.from_string("(((C,B),D),A);"),
        Newick.from_string("(((A,D),B),C);"),
    ]
    expected = Newick.from_string("(((A,B)0,C)2,D);")
    result = main_tree.add_support(bootstraps)
    assert expected == result


def test_newick__add_support__multiple_trees__partially_different_topologies():
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [
        Newick.from_string("(((C,D),A),B);"),
        Newick.from_string("(((A,D),B),C);"),
    ]
    expected = Newick.from_string("(((A,B)1,C)2,D);")
    result = main_tree.add_support(bootstraps)
    assert expected == result


def test_newick__add_support__multiple_trees__two_cladees():
    main_tree = Newick.from_string("((A,B),(C,(D,E)));")
    bootstraps = [
        Newick.from_string("((((C,E),D),A),B);"),
        Newick.from_string("(((A,(C,D)),B),E);"),
    ]
    expected = Newick.from_string("((A,B)1,(C,(D,E)0)1);")
    result = main_tree.add_support(bootstraps)
    assert expected == result


def test_newick__add_support__differing_leaf_names():
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [Newick.from_string("(((C,E),B),A);")]
    with pytest.raises(NewickError):
        main_tree.add_support(bootstraps)


_ADD_SUPPORT_FORMATTING = (
    ("{Support}", "(((A,B)1,C)3,D);"),
    ("{Percentage:.0f}", "(((A,B)33,C)100,D);"),
    ("{Fraction:.2f}", "(((A,B)0.33,C)1.00,D);"),
)


@pytest.mark.parametrize("fmt, expected", _ADD_SUPPORT_FORMATTING)
def test_newick__add_support__formatting(fmt, expected):
    main_tree = Newick.from_string("(((A,B),C),D);")
    bootstraps = [
        Newick.from_string("(((C,D),A),B);"),
        Newick.from_string("(((C,B),A),D);"),
        Newick.from_string("(((A,D),B),C);"),
    ]
    expected = Newick.from_string(expected)
    result = main_tree.add_support(bootstraps, fmt)
    assert expected == result


def test_newick__add_support__unique_names_required():
    main_tree = Newick.from_string("(((A,B),C),A);")
    bootstraps = [Newick.from_string("(((A,B),C),A);")]
    with pytest.raises(NewickError):
        main_tree.add_support(bootstraps)


###############################################################################
###############################################################################
# from_string


def test_newick__parse__minimal_newick__name_only():
    top_node = Newick(name="A")
    assert Newick.from_string("A;") == top_node


def test_newick__parse__single_taxa():
    child_node = Newick(name="Ab")
    top_node = Newick(children=[child_node])
    assert Newick.from_string("(Ab);") == top_node


def test_newick__parse__two_taxa():
    child_node_1 = Newick(name="A")
    child_node_2 = Newick(name="Bc")
    top_node = Newick(children=[child_node_1, child_node_2])
    assert Newick.from_string("(A,Bc);") == top_node


def test_newick__parse__three_taxa():
    child_node_1 = Newick(name="A")
    child_node_2 = Newick(name="Bc")
    child_node_3 = Newick(name="DeF")
    top_node = Newick(children=[child_node_1, child_node_2, child_node_3])
    assert Newick.from_string("(A,Bc,DeF);") == top_node


def test_newick__parse__ignore_whitespace():
    assert Newick.from_string("(A,B);") == Newick.from_string("(A, B);")


def test_newick__parse__missing_semicolon():
    with pytest.raises(NewickParseError):
        Newick.from_string("()")


def test_newick__parse__subnode__single_taxa():
    child_node_1 = Newick(name="A")
    child_node_2a = Newick(name="B")
    child_node_2 = Newick(children=[child_node_2a])
    top_node = Newick(children=[child_node_1, child_node_2])
    assert Newick.from_string("(A,(B));") == top_node


def test_newick__parse__subnode__two_taxa():
    child_node_1 = Newick(name="A")
    child_node_2a = Newick(name="B")
    child_node_2b = Newick(name="C")
    child_node_2 = Newick(children=[child_node_2a, child_node_2b])
    top_node = Newick(children=[child_node_1, child_node_2])
    assert Newick.from_string("(A,(B,C));") == top_node


###########################################################################
###########################################################################
# cmp - white-box, just make sure all properties are compared


def test_newick__cmp__identical():
    node_a = Newick(name="A", length=13, children=[Newick(name="B")])
    node_b = Newick(name="A", length=13, children=[Newick(name="B")])
    assert node_a == node_b


def test_newick__cmp__identical_for_empty_string_length():
    node_a = Newick(name="A", length="", children=[Newick(name="B")])
    node_b = Newick(name="A", length=None, children=[Newick(name="B")])
    assert node_a == node_b


def test_newick__cmp__identical_for_empty_string_name():
    node_a = Newick(name="", length=13, children=[Newick(name="B")])
    node_b = Newick(name=None, length=13, children=[Newick(name="B")])
    assert node_a == node_b


_CMP_NOT_IDENTICAL = (
    Newick(name="B", length=13, children=[Newick(name="B")]),
    Newick(name="A", length=14, children=[Newick(name="B")]),
    Newick(name="A", length=13, children=[]),
    Newick(name="A", length=13, children=[Newick(name="C")]),
    Newick(name="B", length=14, children=[Newick(name="C")]),
)


@pytest.mark.parametrize("node_b", _CMP_NOT_IDENTICAL)
def test_newick__cmp__not_identical(node_b):
    node_a = Newick(name="A", length=13, children=[Newick(name="B")])
    assert node_a != node_b


###############################################################################
###############################################################################
# hash - white-box, just make sure all properties are used


def test_newick__hash__identical():
    node_a = Newick(name="A", length=13, children=[Newick(name="B")])
    node_b = Newick(name="A", length=13, children=[Newick(name="B")])
    assert hash(node_a) == hash(node_b)


_HASH_NOT_IDENTICAL = (
    Newick(name="B", length=13, children=[Newick(name="B")]),
    Newick(name="A", length=14, children=[Newick(name="B")]),
    Newick(name="A", length=13, children=[]),
    Newick(name="A", length=13, children=[Newick(name="C")]),
    Newick(name="B", length=14, children=[Newick(name="C")]),
)


@pytest.mark.parametrize("node_b", _CMP_NOT_IDENTICAL)
def test_newick__hash__not_identical(node_b):
    node_a = Newick(name="A", length=13, children=[Newick(name="B")])
    assert hash(node_a) != hash(node_b)


def test_newick__hash__hashable():
    key_a = Newick(name="A", length=13.7, children=[Newick(name="F")])
    key_b = Newick(name="A", length=13.7, children=[Newick(name="F")])
    assert key_b in {key_a: True}


###############################################################################
###############################################################################
# Malformed newick strings


def test_newick__malformed__unbalanced_parantheses():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C);")


def test_newick__malformed__mismatched_parantheses():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C();")


def test_newick__malformed__missing_parantheses():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C))")


def test_newick__malformed__missing_length():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A:,(B,C));")
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B:,C));")
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C:));")
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C):);")
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,(B,C)):;")


def test_newick__malformed__multiple_lengths():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A:1:2);")


###############################################################################
###############################################################################
# Implicit leafs are not supported (due to problems with ambiguiety)


def test_newick__parse__first_taxa_unnamed():
    with pytest.raises(NewickError):
        Newick.from_string("(,A);")


def test_newick__parse__second_taxa_unnamed():
    with pytest.raises(NewickError):
        Newick.from_string("(A,);")


def test_newick__parse__two_taxa_unnamed():
    with pytest.raises(NewickError):
        Newick.from_string("(,);")


def test_newick__parse__three_taxa_unnamed():
    with pytest.raises(NewickError):
        Newick.from_string("(,,);")


###############################################################################
###############################################################################
# Empty non-leaf nodes are not allowed, as their interpretation is unclear


def test_newick__parse__minimal_newick__implicit_nodes():
    with pytest.raises(NewickParseError):
        Newick.from_string("();")


def test_newick__parse__subnode__empty():
    with pytest.raises(NewickParseError):
        Newick.from_string("(A,());")


###############################################################################
###############################################################################
# The following tests are derived from the wikipedia description of the
# newick format: http://en.wikipedia.org/wiki/Newick_format#Examples


def test_newick__wikipedia_example_1():
    # no nodes are named, this format is not supported here!
    with pytest.raises(NewickError):
        Newick.from_string("(,,(,));")


def test_newick__wikipedia_example_2():
    # leaf nodes are named
    taxa_d = Newick(name="D")
    taxa_c = Newick(name="C")
    taxa_sub = Newick(children=[taxa_c, taxa_d])
    taxa_b = Newick(name="B")
    taxa_a = Newick(name="A")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub])
    assert Newick.from_string("(A,B,(C,D));") == top_node


def test_newick__wikipedia_example_3():
    # all nodes are named
    taxa_d = Newick(name="D")
    taxa_c = Newick(name="C")
    taxa_sub = Newick(children=[taxa_c, taxa_d], name="E")
    taxa_b = Newick(name="B")
    taxa_a = Newick(name="A")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub], name="F")
    assert Newick.from_string("(A,B,(C,D)E)F;") == top_node


def test_newick__wikipedia_example_4():
    # all but root node have a distance to parent
    taxa_d = Newick(length="0.4")
    taxa_c = Newick(length="0.3")
    taxa_sub = Newick(children=[taxa_c, taxa_d], length="0.5")
    taxa_b = Newick(length="0.2")
    taxa_a = Newick(length="0.1")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub])
    assert Newick.from_string("(:0.1,:0.2,(:0.3,:0.4):0.5);") == top_node


def test_newick__wikipedia_example_5():
    # all have a distance to parent
    taxa_d = Newick(length="0.4")
    taxa_c = Newick(length="0.3")
    taxa_sub = Newick(children=[taxa_c, taxa_d], length="0.5")
    taxa_b = Newick(length="0.2")
    taxa_a = Newick(length="0.1")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub], length="0.0")
    assert Newick.from_string("(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;") == top_node


def test_newick__wikipedia_example_6():
    # distances and leaf names (popular)
    taxa_d = Newick(length="0.4", name="D")
    taxa_c = Newick(length="0.3", name="C")
    taxa_sub = Newick(children=[taxa_c, taxa_d], length="0.5")
    taxa_b = Newick(length="0.2", name="B")
    taxa_a = Newick(length="0.1", name="A")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub])
    assert Newick.from_string("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);") == top_node


def test_newick__wikipedia_example_7():
    # distances and all names
    taxa_d = Newick(length="0.4", name="D")
    taxa_c = Newick(length="0.3", name="C")
    taxa_sub = Newick(children=[taxa_c, taxa_d], length="0.5", name="E")
    taxa_b = Newick(length="0.2", name="B")
    taxa_a = Newick(length="0.1", name="A")
    top_node = Newick(children=[taxa_a, taxa_b, taxa_sub], name="F")
    assert Newick.from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;") == top_node


def test_newick__wikipedia_example_8():
    # a tree rooted on a leaf node (rare)
    taxa_b = Newick(length="0.2", name="B")
    taxa_c = Newick(length="0.3", name="C")
    taxa_d = Newick(length="0.4", name="D")
    node_e = Newick(length="0.5", name="E", children=[taxa_c, taxa_d])
    node_f = Newick(length="0.1", name="F", children=[taxa_b, node_e])
    node_a = Newick(name="A", children=[node_f])
    assert Newick.from_string("((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;") == node_a


###############################################################################
###############################################################################
# str / repr


def test_newick__str__non_string_name():
    node = Newick(children=[Newick(name=17, length="1.3")])
    assert str(node) == "(17:1.3);"


def test_newick__str__non_string_length():
    node = Newick(children=[Newick(name="Foo", length=1.3)])
    assert str(node) == "(Foo:1.3);"


def test_newick__str__repr_equal_to_str():
    node_a = Newick(name="A", length="123")
    node_b = Newick(name="B", length="4567")
    top_node = Newick(children=[node_a, node_b])
    assert str(top_node) == "(A:123,B:4567);"


def test_newick__str__single_leaf_should_not_be_followed_by_comma():
    node = Newick(name="A")
    top_node = Newick(children=[node])
    assert str(top_node) == "(A);"


_STR_EQUALITY = (
    # 2. leaf nodes are named
    "(A,B,(C,D));",
    # 3. all nodes are named
    "(A,B,(C,D)E)F;",
    # 4. all but root node have a distance to parent
    "(:0.1,:0.2,(:0.3,:0.4):0.5);",
    # 5. all have a distance to parent
    "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
    # 6. distances and leaf names (popular)
    "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
    # 7. distances and all names
    "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
)


@pytest.mark.parametrize("nwk_str", _STR_EQUALITY)
def test_newick__wikipedia_examples__str_equality(nwk_str):
    nodes = Newick.from_string(nwk_str)
    result = str(nodes)
    assert result == nwk_str


###############################################################################
###############################################################################
# Immutability

_PROPERTIES_ARE_IMMUTABLE = (
    ("name", "foo"),
    ("length", "13"),
    ("children", []),
    ("foobar", True),
)


@pytest.mark.parametrize("name, value", _PROPERTIES_ARE_IMMUTABLE)
def test_newick__properties_are_immutable(name, value):
    node = Newick(name="A", length=3, children=[Newick(name="B")])
    with pytest.raises(NotImplementedError):
        setattr(node, name, value)


@pytest.mark.parametrize("name", ("name", "length", "children", "foobar"))
def test_newick__properties_cannot_be_deleted(name):
    node = Newick(name="A", length=3, children=[Newick(name="B")])
    with pytest.raises(NotImplementedError):
        delattr(node, name)
