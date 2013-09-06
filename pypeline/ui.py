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
import datetime
import collections

import pypeline.logger
from pypeline.node import MetaNode
from pypeline.common.console import \
     print_msg, \
     print_debug, \
     print_info, \
     print_err, \
     print_warn, \
     print_disabled


def print_node_tree(graph, collapse = True, verbose = True):
    print_msg(datetime.datetime.now().strftime("%F %T"))
    print_msg("Pipeline%s" % _describe_nodes(graph, graph.iterflat()))

    logfile = pypeline.logger.get_logfile()
    if logfile:
        print_debug("  Log-file located at %r" % (logfile,))

    if verbose:
        _print_sub_nodes(graph, graph, collapse, "  ")
    else:
        _print_running_nodes(graph)


def _print_running_nodes(graph):
    running = []
    for node in graph.iterflat():
        if graph.get_node_state(node) is graph.RUNNING:
            running.append(node)

    for node in sorted(running, key = str):
        print_info("  - %s" % node)
    print_info()


def _print_sub_nodes(graph, nodes, collapse, prefix = ""):
    viable_nodes, dead_nodes = [], []
    for node in nodes:
        if collapse and (graph.get_node_state(node) == graph.DONE):
            dead_nodes.append(node)
        else:
            viable_nodes.append(node)
    viable_nodes.sort(key = str)

    for node in viable_nodes:
        description = "%s%s %s" % (prefix, _get_runable_prefix(graph, node), node)
        if node.subnodes:
            description += _describe_nodes(graph, node.subnodes)

        print_func = _get_print_function(graph, node)
        print_func(description)

        is_last_node = (node == viable_nodes[-1]) and not dead_nodes
        current_prefix = prefix + ("  " if is_last_node else "| ")

        if node.dependencies:
            if collapse and _collapse_node(graph, node.dependencies):
                description = "+ %i dependencies hidden ..." \
                    % _count_dependencies(node.dependencies | node.subnodes)

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
    def _do_count_dependencies(dependencies, observed):
        novel_nodes = (dependencies - observed)
        observed.update(dependencies)
        for node in novel_nodes:
            _do_count_dependencies(node.dependencies | node.subnodes, observed)
        return observed

    count = 0
    for node in _do_count_dependencies(set(dependencies), set()):
        if not isinstance(node, MetaNode):
            count += 1
    return count


def _describe_nodes(graph, nodes):
    states = collections.defaultdict(int)
    def count_states(node_lst, recurse = True):
        for node in node_lst:
            if not isinstance(node, MetaNode):
                states[graph.get_node_state(node)] += 1
            elif recurse:
                count_states(node.subnodes, recurse = False)
    count_states(nodes)

    fields = [("running",  states[graph.RUNNING]),
              ("outdated", states[graph.OUTDATED]),
              ("failed",   states[graph.ERROR])]

    line = [""]
    for (name, value) in fields:
        if value:
            line.append("%i %s" % (value, name))

    line.append("%i done of %i tasks" \
                    % (states[graph.DONE],
                       sum(states.values())))

    return ", ".join(line)


def _collapse_node(graph, dependencies):
    """Returns true if a node may be collapsed in the dependency graph."""
    if all((graph.get_node_state(node) == graph.DONE) for node in dependencies):
        return (_count_dependencies(dependencies) > 2)

    return False


def _get_print_function(graph, node):
    for subnode in node.subnodes:
        if graph.get_node_state(subnode) == graph.RUNNING:
            return print_info

    state = graph.get_node_state(node)
    if state is graph.RUNNING:
        return print_info
    elif state is graph.DONE:
        return print_disabled
    elif state is graph.OUTDATED:
        return print_warn
    elif state is graph.ERROR:
        return print_err
    else:
        return print_msg


def _get_runable_prefix(graph, node):
    """Returns either 'R' or '+', dependening on the state of the node. If the node, or any
    of its subnodes, are runable, then 'R' is returned, otherwise '+' is returned. This is
    used to decorate the dependency graph."""
    if graph.get_node_state(node) in (graph.RUNNING, graph.RUNABLE):
        return "R"

    for subnode in node.subnodes:
        if graph.get_node_state(subnode) in (graph.RUNNING, graph.RUNABLE):
            return "R"

    return "+"
