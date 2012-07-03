"""Functions relating to the CLI interface."""
from __future__ import print_function

from node import MetaNode


def _do_print_color(*vargs, **kwargs):
    """Utility function: Prints using shell colors."""
    colorcode = kwargs.pop("colorcode")
    text = ["\033[00;%im%s\033[00m" % (colorcode, arg) for arg in vargs]

    print(*text, **kwargs)


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



def print_node_tree(top_nodes, running):
    print_msg("Pipeline (%i nodes running):" % (len(running),))
    _print_sub_nodes(top_nodes, running, "   ")
        

def _print_sub_nodes(nodes, running, prefix = ""):
    nodes = list(nodes)
    nodes.sort(key = str)

    for node in nodes:
        if node.is_runable and not node.is_done or node in running:
            description = prefix + "R " + str(node)
        else:
            description = prefix + "+ " + str(node)

        if isinstance(node, MetaNode):
            active, done, total = 0, 0, 0
            for subnode in node.subnodes:
                total += 1
                if subnode in running:
                    active += 1
                elif subnode.is_done:
                    done += 1
        
            description += " %i running, %i done of %i subnodes" % (active, done, total)

        print_func = _get_print_function(node, running)
        print_func(description)
        current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

        dependencies = _collect_dependencies(node)
        if dependencies:
            _print_sub_nodes(dependencies, running, current_prefix + "   ")
        else:
            print_func(current_prefix)


def _get_print_function(node, running):
    if isinstance(node, MetaNode):
        for subnode in node.subnodes:
            if subnode in running:
                return print_info

    if node in running:
        return print_info
    elif node.is_done:
        return print_disabled
    elif node.is_outdated:
        return print_warn
    else:
        return print_msg


def _collect_dependencies(node):
    """For a regular node, this function returns the subnodes, which are taken
    to be the dependencies of that node. For a MetaNode, the subnodes are 
    considered part of that Node, and hence the dependencies of _those_ nodes
    are returned."""
    if not isinstance(node, MetaNode):
        return node.subnodes

    dependencies = set()
    for subnode in node.subnodes:
        if isinstance(subnode, MetaNode):
            dependencies.add(subnode)
        else:
            dependencies.update(subnode.subnodes)

    return dependencies
