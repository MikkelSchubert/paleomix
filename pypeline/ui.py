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
import time
import datetime
import optparse

import pypeline.nodegraph
import pypeline.logger
from pypeline.node import MetaNode
from pypeline.common.console import \
     print_msg, \
     print_debug, \
     print_info, \
     print_err, \
     print_warn, \
     print_disabled
from pypeline.common.utilities import \
     group_by_pred



def add_optiongroup(parser, default = "quiet"):
    """Adds an option-group to an OptionParser object, with options
    pertaining to logging. Note that 'initialize' expects the config
    object to have these options."""
    group = optparse.OptionGroup(parser, "Progress reporting")
    group.add_option("--progress-ui", default = default, type = "choice",
                     choices = ("verbose", "quiet", "progress"),
                     help = "Select method for displaying the progress of the pipeline: "
                            "'verbose' = Full dependency tree at every change; "
                            "'quiet' = Display only currently running nodes; "
                            "'progress' = Display changes in state. "
                            "[Default is '%default']")
    parser.add_option_group(group)


def get_ui(ui_name):
    """Returns a UI instance by name, using the choices allowed by
    the 'add_optiongroup' function. See keys in 'UI_TYPES'."""
    ui_name = ui_name.title()
    if ui_name not in UI_TYPES:
        raise ValueError("Unknown UI type %r" % (ui_name,))
    return UI_TYPES[ui_name]()




class BaseUI:
    """UI base class.

    Can be initialized, but does nothing but collect stats about
    the pipeline. Subclasses should override at least one of
    (but still call the BaseUI function) the functions 'flush',
    'finalize', and/or 'state_changed'.

    In addition, the class contains the following properties:
      - states  -- List containing the observed number of states
                   for a state-value corresponding to the index
      - threads -- Est. number of threads used by running nodes.

    These properties should be treated as read-only.
    """

    def __init__(self):
        """Basic initializer; must be called in subclasses."""
        self.states  = []
        self.threads = 0


    def flush(self):
        """Called by the user of the UI to ensure that the UI to print
        the current state of the pipeline / changes to pipeline / etc."""
        pass


    def finalize(self):
        """Called by the pipeline at the termination of a run. By default,
        this function prints the location of the log-file if one was created
        during the run (e.g. if there were errors)."""
        logfile = pypeline.logger.get_logfile()
        if logfile:
            print_debug("Log-file located at %r" % (logfile,))

        if self.states[self.ERROR]:
            print_err("Done; but errors were detected ...")
        else:
            print_info("Done ...")


    def refresh(self, nodegraph):
        """Called when the nodegraph has refreshed, causing state-counts
        to be recalculated."""
        self.states, self.threads \
          = self._count_states(nodegraph, nodegraph.iterflat())


    def state_changed(self, node, old_state, new_state, _is_primary):
        """Observer function for NodeGraph; counts states for non-meta nodes."""
        if not isinstance(node, MetaNode):
            self.states[old_state] -= 1
            self.states[new_state] += 1
            if old_state == self.RUNNING:
                self.threads -= node.threads
            elif new_state == self.RUNNING:
                self.threads += node.threads


    @classmethod
    def _count_states(self, nodegraph, nodes, meta = False):
        """Counts the number of each state observed for a set of nodes, and
        returns these as a list, as well as the estimated number of threads
        being used by running nodes.

        If 'meta' is true, these are considered to be the subnodes of a
        MetaNode; in that case, MetaNode(s) themselves are not counted, but
        their subnodes are counted for the first level encountered."""
        states  = [0] * nodegraph.NUMBER_OF_STATES
        threads = []

        def inc_states(c_nodes, depth):
            for node in c_nodes:
                is_meta = isinstance(node, MetaNode)

                if depth and is_meta:
                    inc_states(node.subnodes, depth - 1)
                elif meta or not is_meta:
                    state = nodegraph.get_node_state(node)
                    states[state] += 1
                    if state == nodegraph.RUNNING:
                        threads.append(node.threads)

            return states

        return inc_states(nodes, (1 if meta else 0)), sum(threads)


    @classmethod
    def _describe_states(cls, states, threads = 0):
        """Returns a human readable summary of the states
        given to the function. 'states' is expected to be
        a list/tuple such as that produced by the
        '_count_states' function."""
        run_tmpl = "running"
        if threads > 1:
            run_tmpl = "%s using ~%i threads" % (run_tmpl, threads)
        elif threads == 1:
            run_tmpl = "%s using ~1 thread" % (run_tmpl,)

        fields = [(run_tmpl,   states[cls.RUNNING]),
                  ("outdated", states[cls.OUTDATED]),
                  ("failed",   states[cls.ERROR])]

        line = []
        for (name, value) in fields:
            if value:
                line.append("%i %s" % (value, name))

        line.append("%i done of %i tasks" \
                    % (states[cls.DONE],
                       sum(states)))

        return ", ".join(line)


    DONE     = pypeline.nodegraph.NodeGraph.DONE
    RUNNING  = pypeline.nodegraph.NodeGraph.RUNNING
    RUNABLE  = pypeline.nodegraph.NodeGraph.RUNABLE
    QUEUED   = pypeline.nodegraph.NodeGraph.QUEUED
    OUTDATED = pypeline.nodegraph.NodeGraph.OUTDATED
    ERROR    = pypeline.nodegraph.NodeGraph.ERROR




class VerboseUI(BaseUI):
    def __init__(self):
        BaseUI.__init__(self)
        self._graph = None


    def flush(self):
        """See BaseUI.flush."""
        BaseUI.flush(self)
        self._print_header(self.states, self.threads)
        self._print_sub_nodes(self._graph, self._graph)


    def refresh(self, nodegraph):
        """See BaseUI.refresh."""
        BaseUI.refresh(self, nodegraph)
        self._graph = nodegraph


    @classmethod
    def _print_header(cls, states, threads):
        print_msg(datetime.datetime.now().strftime("%F %T"))
        print_msg("Pipeline; %s:" % cls._describe_states(states, threads))
        logfile = pypeline.logger.get_logfile()
        if logfile:
            print_debug("  Log-file located at %r" % (logfile,))


    @classmethod
    def _print_sub_nodes(cls, nodegraph, nodes, prefix = "  "):
        viable_nodes, dead_nodes = \
          group_by_pred(lambda node: (nodegraph.get_node_state(node) != cls.DONE), nodes)
        viable_nodes.sort(key = str)

        for node in viable_nodes:
            runable     = cls._get_runable_prefix(nodegraph, node)
            description = "%s%s %s" % (prefix, runable, node)
            if node.subnodes:
                states, threads = cls._count_states(nodegraph, node.subnodes, True)
                description = "%s (%s)" % (description, cls._describe_states(states, threads))

            print_func = cls._get_print_function(nodegraph, node)
            print_func(description)

            is_last_node = (node == viable_nodes[-1]) and not dead_nodes
            current_prefix = prefix + ("  " if is_last_node else "| ")

            if node.dependencies:
                if cls._collapse_node(nodegraph, node.dependencies):
                    description = "+ %i dependencies hidden ..." \
                      % cls._count_dependencies(node.dependencies | node.subnodes)

                    print_disabled(current_prefix + description)
                    print_disabled(current_prefix)
                else:
                    cls._print_sub_nodes(nodegraph, node.dependencies, current_prefix + "  ")
            else:
                print_func(current_prefix)

        if dead_nodes:
            print_disabled(prefix + "+ %i dependencies hidden ..." \
                           % cls._count_dependencies(dead_nodes))
            print_disabled(prefix)


    @classmethod
    def _get_runable_prefix(cls, nodegraph, node):
        """Returns either 'R' or '+', dependening on the state of the node. If the node, or any
        of its subnodes, are runable, then 'R' is returned, otherwise '+' is returned. This is
        used to decorate the dependency graph."""
        if nodegraph.get_node_state(node) in (cls.RUNNING, cls.RUNABLE):
            return "R"

        for subnode in node.subnodes:
            if nodegraph.get_node_state(subnode) in (cls.RUNNING, cls.RUNABLE):
                return "R"

        return "+"


    @classmethod
    def _count_dependencies(cls, dependencies):
        """Recursively counts dependencies / subnodes in a subtree,
        excluding MetaNodes from the total. Node that are depended
        upon by multiple other nodes are only counted once."""
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


    @classmethod
    def _collapse_node(cls, graph, dependencies):
        """Returns true if a node may be collapsed in the dependency graph."""
        if all((graph.get_node_state(node) == graph.DONE) for node in dependencies):
            return (cls._count_dependencies(dependencies) > 2)

        return False


    @classmethod
    def _get_print_function(cls, graph, node):
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
        return print_msg




class QuietUI(VerboseUI):
    """A more quiet progress UI, relative to the Verbose UI:
    Prints a summary, and the list of running nodes every
    time 'flush' is called."""

    def __init__(self):
        VerboseUI.__init__(self)
        self._running_nodes = []


    def flush(self):
        """See BaseUI.flush."""
        if not self._running_nodes:
            VerboseUI.flush(self)
            return

        BaseUI.flush(self)
        self._print_header(self.states, self.threads)
        for node in sorted(self._running_nodes, key = str):
            print_info("  - %s" % node)
        print_info()


    def state_changed(self, node, old_state, new_state, is_primary):
        """See BaseUI.state_changed."""
        BaseUI.state_changed(self, node, old_state, new_state, is_primary)
        if not is_primary:
            # Avoid printing updates for nodes that changed because another
            # node changed. E.g. don't print all the nodes that FAILED
            # because a dependency FAILED.
            return

        if old_state == self.RUNNING:
            self._running_nodes.remove(node)
        elif new_state == self.RUNNING:
            self._running_nodes.append(node)



class ProgressUI(BaseUI):
    """Progress based UI: Prints nodes when they start running, when they finish
    running, or when they fail running. Changes to state resulting from the
    above is not printed. Every 25th update is followed by a summary of the
    current total progress."""

    # Print a summery of the current state very N events
    _SUMMARY_EVERY = 25

    def __init__(self):
        self._refresh_count = ProgressUI._SUMMARY_EVERY
        self._runtimes = {}
        BaseUI.__init__(self)


    def refresh(self, nodegraph):
        """See BaseUI.refresh."""
        BaseUI.refresh(self, nodegraph)
        self._print_summary()


    def state_changed(self, node, old_state, new_state, is_primary):
        """See BaseUI.state_changed."""
        BaseUI.state_changed(self, node, old_state, new_state, is_primary)
        if is_primary and (new_state in self._DESCRIPTIONS):
            if new_state in (self.RUNNING, self.DONE, self.ERROR):
                self._print_state(node, new_state)

            self._refresh_count -= 1
            if (self._refresh_count <= 0) or (new_state == self.ERROR):
                self._refresh_count = ProgressUI._SUMMARY_EVERY
                self._print_summary()


    def _print_summary(self):
        """Prints a summary of the pipeline progress."""
        time_label  = datetime.datetime.now().strftime("%T")
        description = self._describe_states(self.states, self.threads)
        print_msg("\n%s Pipeline: %s" % (time_label, description))
        logfile = pypeline.logger.get_logfile()
        if logfile:
            print_debug("Log-file located at %r" % (logfile,))


    def _print_state(self, node, new_state):
        state_label, print_func = self._DESCRIPTIONS[new_state]
        if new_state == self.RUNNING:
            self._runtimes[node] = time.time()
        elif new_state in (self.RUNNING, self.DONE, self.ERROR):
            state_label = "%s (%s)" % (state_label, self._get_runtime(node))

        time_label = datetime.datetime.now().strftime("%T")
        print_func("%s %s: %s" % (time_label, state_label, node))


    def _get_runtime(self, node):
        current_time = time.time()
        runtime = int(current_time - self._runtimes.pop(node, current_time))
        if runtime >= 3600:
            fmt = "{hours}:{mins:02}:{secs:02}s"
        elif runtime >= 60:
            fmt = "{mins}:{secs:02}s"
        else:
            fmt = "{secs}s"

        return fmt.format(hours = runtime // 3600,
                          mins  = (runtime // 60) % 60,
                          secs  = (runtime % 60))


    _DESCRIPTIONS = {
        BaseUI.DONE    : ("Finished", print_disabled),
        BaseUI.RUNNING : ("Started",  print_info),
        BaseUI.ERROR   : ("Failed",   print_err),
    }


# Different types of UIs
UI_TYPES = {
    "Verbose"  : VerboseUI,
    "Quiet"    : QuietUI,
    "Progress" : ProgressUI,
    }
