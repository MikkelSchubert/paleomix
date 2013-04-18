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
"""Functions relating to the CLI interface."""
from __future__ import print_function

import sys
import datetime
import collections

from pypeline.nodegraph import NodeGraph
from pypeline.node import MetaNode


def _do_print_color(*vargs, **kwargs):
    """Utility function: Prints using shell colors."""
    colorcode = kwargs.pop("colorcode")
    destination = kwargs.pop("file", sys.stdout)

    # No colors if output is redirected (e.g. less).
    if destination.isatty():
        vargs = ["\033[00;%im%s\033[00m" % (colorcode, arg) for arg in vargs]

    print(*vargs, file = destination, **kwargs)


def print_msg(*vargs, **kwargs):
    """Equivalent to print. Currently does not apply a color to the text"""
    print(*vargs, **kwargs)


def print_info(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (green)."""
    _do_print_color(*vargs, colorcode = 32, **kwargs)


def print_err(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (red)."""
    _do_print_color(*vargs, colorcode = 31, **kwargs)


def print_warn(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (yellow)."""
    _do_print_color(*vargs, colorcode = 33, **kwargs)


def print_disabled(*vargs, **kwargs):
    """Equivalent to print, but prints using shell colorcodes (gray)."""
    _do_print_color(*vargs, colorcode = 30, **kwargs)



def print_node_tree(graph, collapse = True, verbose = True):
    print_msg(datetime.datetime.now().strftime("%F %T"))
    print_msg("Pipeline%s" % _describe_nodes(graph, graph.iterflat()))

    if verbose:
        _print_sub_nodes(graph, graph, collapse, "  ")
    else:
        _print_running_nodes(graph)


def _print_running_nodes(graph):
    running = []
    for node in graph.iterflat():
        if graph.get_node_state(node) is NodeGraph.RUNNING:
            running.append(node)

    for node in sorted(running, key = str):
        print_info("  - %s" % node)
    print_info()


def _print_sub_nodes(graph, nodes, collapse, prefix = ""):
    viable_nodes, dead_nodes = [], []
    for node in nodes:
        if collapse and (graph.get_node_state(node) == NodeGraph.DONE):
            dead_nodes.append(node)
        else:
            viable_nodes.append(node)
    viable_nodes.sort(key = str)

    for node in viable_nodes:
        description = "%s%s %s" % (prefix, _get_runable_prefix(graph, node), node)
        if node.subnodes:
            description += _describe_nodes(graph, node.subnodes, count_meta = True)

        print_func = _get_print_function(graph, node)
        print_func(description)

        is_last_node = (node == viable_nodes[-1]) and not dead_nodes
        current_prefix = prefix + ("  " if is_last_node else "| ")

        if node.dependencies:
            if collapse and _collapse_node(graph, node.dependencies):
                description = "+ %i dependencies hidden ..." \
                    % _count_dependencies(node.dependencies)

                print_disabled(current_prefix + description)
                print_disabled(current_prefix)
            else:
                _print_sub_nodes(graph, node.dependencies, collapse, current_prefix + "  ")
        else:
            print_func(current_prefix)

    if dead_nodes:
        print_disabled(prefix + "+ %i dependencies hidden ..." \
                           % _count_dependencies(dead_nodes))
        print_disabled(prefix)


def _count_dependencies(dependencies):
    counter = len(dependencies)
    for node in dependencies:
        counter += _count_dependencies(node.dependencies)

    return counter


def _describe_nodes(graph, nodes, count_meta = False):
    states = collections.defaultdict(int)
    for node in nodes:
        if count_meta or not isinstance(node, MetaNode):
            states[graph.get_node_state(node)] += 1

    fields = [("running",  states[NodeGraph.RUNNING]),
              ("outdated", states[NodeGraph.OUTDATED]),
              ("failed",   states[NodeGraph.ERROR])]

    line = [""]
    for (name, value) in fields:
        if value:
            line.append("%i %s" % (value, name))

    line.append("%i done of %i tasks" \
                    % (states[NodeGraph.DONE], 
                       sum(states.values())))

    return ", ".join(line)
    

def _collapse_node(graph, dependencies):
    """Returns true if a node may be collapsed in the dependency graph."""
    if all((graph.get_node_state(node) == NodeGraph.DONE) for node in dependencies):
        return (_count_dependencies(dependencies) > 2)

    return False


def _get_print_function(graph, node):
    for subnode in node.subnodes:
        if graph.get_node_state(subnode) == NodeGraph.RUNNING:
            return print_info

    state = graph.get_node_state(node)
    if state is NodeGraph.RUNNING:
        return print_info
    elif state is NodeGraph.DONE:
        return print_disabled
    elif state is NodeGraph.OUTDATED:
        return print_warn
    elif state is NodeGraph.ERROR:
        return print_err
    else:
        return print_msg


def _get_runable_prefix(graph, node):
    """Returns either 'R' or '+', dependening on the state of the node. If the node, or any
    of its subnodes, are runable, then 'R' is returned, otherwise '+' is returned. This is 
    used to decorate the dependency graph."""
    if graph.get_node_state(node) in (NodeGraph.RUNNING, NodeGraph.RUNABLE):
        return "R"
    
    for subnode in node.subnodes:
        if graph.get_node_state(subnode) in (NodeGraph.RUNNING, NodeGraph.RUNABLE):
            return "R"
    
    return "+"
