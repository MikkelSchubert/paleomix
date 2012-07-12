"""Functions relating to the CLI interface."""
from __future__ import print_function

import sys

from node import MetaNode


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
    total, done, running, failed = 0, 0, 0, 0
    for node in graph.iterflat():
        total += 1
        if node.state == node.DONE:
            done += 1
        elif node.state == node.RUNNING:
            running += 1
        elif node.state == node.ERROR:
            failed  += 1
            

    print_msg("Pipeline, \t%i running, %i done, %i failed of %i nodes:" \
                  % (running, done, failed, total))
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
        if node.state in (node.RUNNING, node.RUNABLE):
            description = prefix + "R " + str(node)
        else:
            description = prefix + "+ " + str(node)

        # FIXME
        if isinstance(node.task, MetaNode):
            active, done, outdated, total = 0, 0, 0, 0
            for subnode in node.subnodes:
                total += 1
                if subnode.state == subnode.RUNNING:
                    active += 1
                elif subnode.state == subnode.DONE:
                    done += 1
                elif subnode.state == subnode.OUTDATED:
                    outdated += 1
         
            description += " %i running, %i outdated, %i done of %i subnodes" \
                % (active, outdated, done, total)

        print_func = _get_print_function(node)
        print_func(description)

        is_last_node = (node == viable_nodes[-1]) and not dead_nodes
        current_prefix = prefix + ("  " if is_last_node else "|  ")

        dependencies = _collect_dependencies(node)
        if dependencies:
            if collapse and _collapse_node(dependencies):
                description = "+ %i nodes hidden ..." % _count_subnodes(dependencies)
                print_disabled(current_prefix + description)
                print_disabled(current_prefix)
            else:
                _print_sub_nodes(dependencies, collapse, current_prefix + "   ")
        else:
            print_func(current_prefix)

    if dead_nodes:
        print_disabled(prefix + "+ %i nodes hidden ..." % _count_subnodes(dead_nodes))
        print_disabled(prefix)


def _count_subnodes(nodes):
    counter = len(nodes)
    for node in nodes:
        counter += _count_subnodes(_collect_dependencies(node))
    return counter


def _collapse_node(dependencies):
    for subnode in dependencies:
        if subnode.state != subnode.DONE:
            return False

    return True


def _get_print_function(node):
    # FIXME
    if isinstance(node.task, MetaNode):
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


def _collect_dependencies(node):
    """For a regular node, this function returns the subnodes, which are taken
    to be the dependencies of that node. For a MetaNode, the subnodes are 
    considered part of that Node, and hence the dependencies of _those_ nodes
    are returned."""
    # FIXME
    if not isinstance(node.task, MetaNode):
        return node.subnodes

    dependencies = set()
    for subnode in node.subnodes:
        for dependency in subnode.subnodes:
            if dependency not in node.subnodes:
                dependencies.add(dependency)
                
    return dependencies
