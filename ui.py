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



def print_node_tree(nodes, states, collapse = True):
    print_msg("Pipeline (%i nodes running):" % (states.running_tasks(),))
    _print_sub_nodes(nodes, states, collapse, "   ")
        

def _print_sub_nodes(nodes, states, collapse, prefix = ""):
    nodes = list(nodes)
    nodes.sort(key = str)

    for node in nodes:
        state = states.get_state(node)
        if state in (states.RUNNING, states.RUNABLE):
            description = prefix + "R " + str(node)
        else:
            description = prefix + "+ " + str(node)

        if isinstance(node, MetaNode):
            active, done, total = 0, 0, 0
            for subnode in node.subnodes:
                total += 1
                substate = states.get_state(subnode)
                if substate == states.RUNNING:
                    active += 1
                elif substate == states.DONE:
                    done += 1
        
            description += " %i running, %i done of %i subnodes" % (active, done, total)

        print_func = _get_print_function(node, states)
        print_func(description)
        current_prefix = prefix + ("  " if (node == nodes[-1]) else "|  ")

        dependencies = _collect_dependencies(node)
        if dependencies:
            if collapse and _collapse_node(dependencies, states):
                print_disabled(current_prefix + "+ ...")
                print_disabled(current_prefix)
            else:
                _print_sub_nodes(dependencies, states, collapse, current_prefix + "   ")
        else:
            print_func(current_prefix)


def _collapse_node(dependencies, states):
    for subnode in dependencies:
        if states.get_state(subnode) != states.DONE:
            return False

    return True


def _get_print_function(node, states):
    if isinstance(node, MetaNode):
        for subnode in node.subnodes:
            if states.get_state(subnode) == states.RUNNING:
                return print_info

    state = states.get_state(node)
    if state is states.RUNNING:
        return print_info
    elif state is states.DONE:
        return print_disabled
    elif state is states.OUTDATED:
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
        for dependency in subnode.subnodes:
            if dependency not in node.subnodes:
                dependencies.add(dependency)
                
    return dependencies
