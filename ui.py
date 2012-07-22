"""Functions relating to the CLI interface."""
from __future__ import print_function

import sys
import collections

from pypeline.taskgraph import TaskGraph


def _do_print_color(*vargs, **kwargs):
    """Utility function: Prints using shell colors."""
    colorcode = kwargs.pop("colorcode")
    # No colors if stdout is redirected (e.g. less). 
    if sys.stdout.isatty():
        vargs = ["\033[00;%im%s\033[00m" % (colorcode, arg) for arg in vargs]

    print(*vargs, **kwargs)


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



def print_node_tree(graph, collapse = True):
    print_msg("Pipeline,%s" % _describe_nodes(graph.iterflat()))
    _print_sub_nodes(graph, collapse, "   ")
        

def _print_sub_nodes(nodes, collapse, prefix = ""):
    viable_nodes, dead_nodes = [], []
    for node in nodes:
        if collapse and (node.state == node.DONE):
            dead_nodes.append(node)
        else:
            viable_nodes.append(node)
    viable_nodes.sort(key = str)

    for node in viable_nodes:
        description = "%s%s %s" % (prefix, _get_runable_prefix(node), node)
        if node.subnodes:
            description += _describe_nodes(node.subnodes)
            
        print_func = _get_print_function(node)
        print_func(description)

        is_last_node = (node == viable_nodes[-1]) and not dead_nodes
        current_prefix = prefix + ("   " if is_last_node else "|  ")

        if node.dependencies:
            if collapse and _collapse_node(node.dependencies):
                description = "+ %i dependencies hidden ..." \
                    % _count_dependencies(node.dependencies)

                print_disabled(current_prefix + description)
                print_disabled(current_prefix)
            else:
                _print_sub_nodes(node.dependencies, collapse, current_prefix + "   ")
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


def _describe_nodes(nodes):
    states = collections.defaultdict(int)
    for node in nodes:
        states[node.state] += 1

    return " %(running)i running, %(outdated)i outdated, %(failed)i failed, %(done)i done of %(total)i nodes" \
        % {"running"  : states[TaskGraph.Node.RUNNING], 
           "outdated" : states[TaskGraph.Node.OUTDATED], 
           "done"     : states[TaskGraph.Node.DONE], 
           "failed"   : states[TaskGraph.Node.ERROR], 
           "total"    : sum(states.values())}


def _collapse_node(dependencies):
    """Returns true if a node may be collapsed in the dependency graph."""
    if all((node.state == node.DONE) for node in dependencies):
        return (_count_dependencies(dependencies) > 2)

    return False


def _get_print_function(node):
    for subnode in node.subnodes:
        if subnode.state == subnode.RUNNING:
            return print_info

    if node.state is node.RUNNING:
        return print_info
    elif node.state is node.DONE:
        return print_disabled
    elif node.state is node.OUTDATED:
        return print_warn
    else:
        return print_msg


def _get_runable_prefix(node):
    """Returns either 'R' or '+', dependening on the state of the node. If the node, or any
    of its subnodes, are runable, then 'R' is returned, otherwise '+' is returned. This is 
    used to decorate the dependency graph."""
    if node.state in (node.RUNNING, node.RUNABLE):
        return "R"
    
    for subnode in node.subnodes:
        if subnode.state in (node.RUNNING, node.RUNABLE):
            return "R"
    
    return "+"
