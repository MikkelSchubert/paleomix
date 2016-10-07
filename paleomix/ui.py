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
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
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
import multiprocessing
import optparse
import os
import select
import sys
import termios
import time
import tty

import paleomix.nodegraph
import paleomix.logger
import paleomix.common.text as text

from paleomix.common.console import \
    print_msg, \
    print_debug, \
    print_info, \
    print_err, \
    print_disabled, \
    print_warn


def add_optiongroup(parser, ui_default="running", color_default="on"):
    """Adds an option-group to an OptionParser object, with options
    pertaining to logging. Note that 'initialize' expects the config
    object to have these options."""
    group = optparse.OptionGroup(parser, "Progress reporting")
    group.add_option("--progress-ui", default=ui_default, type="choice",
                     choices=("running", "progress", "summary",
                              "verbose", "quiet"),
                     help="Select method for displaying the progress of the "
                          "pipeline: 'running' = Display only currently "
                          "running nodes; 'progress' = Display changes in "
                          "state; 'summary'; one-line summary only. "
                          "[Default is '%default']")
    group.add_option("--ui-colors", default=color_default,
                     choices=("on", "off", "force"),
                     help="Enable, disable, or force the use of color codes "
                          "when printing the command-line UI. Unless forced, "
                          "colors will only be printed if STDOUT is a TTY "
                          "[Default is '%default']")
    parser.add_option_group(group)


def get_ui(ui_name):
    """Returns a UI instance by name, using the choices allowed by
    the 'add_optiongroup' function. See keys in 'UI_TYPES'."""
    ui_name = ui_name.title()
    if ui_name not in UI_TYPES:
        raise ValueError("Unknown UI type %r" % (ui_name,))
    return UI_TYPES[ui_name]()


def set_ui_colors(choice):
    import paleomix.common.console as console
    choice = choice.lower()
    if choice == "on":
        console.set_color_output(console.COLORS_ON)
    elif choice == "off":
        console.set_color_output(console.COLORS_OFF)
    elif choice == "force":
        console.set_color_output(console.COLORS_FORCED)
    else:
        raise ValueError("Unknown color setting %r" % (choice,))


class CommandLine(object):
    def __init__(self):
        self._tty_settings = None

    def __enter__(self):
        assert not self._tty_settings
        # False if the pipeline is being piped somewhere
        if sys.stdin.isatty() and sys.stdout.isatty():
            # False if the process is running in the background
            if os.getpgrp() == os.tcgetpgrp(sys.stdout.fileno()):
                try:
                    # Store old settings
                    self._tty_settings = termios.tcgetattr(sys.stdin)
                    # Disable echo
                    tty.setcbreak(sys.stdin.fileno())
                except tty.error:
                    pass  # Silently ignore failures

        return self

    def __exit__(self, _type, _value, _traceback):
        if self._tty_settings:
            # Restore settings (re-enable echo)
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, self._tty_settings)

    def process_key_presses(self, nodegraph, max_threads, ui):
        if not self._tty_settings:
            return max_threads

        help_printed = False
        old_max_threads = max_threads
        while self.poll_stdin():
            character = sys.stdin.read(1)
            if character == "+":
                max_threads = min(multiprocessing.cpu_count(), max_threads + 1)
            elif character == "-":
                max_threads = max(1, max_threads - 1)
            elif character in "lL":
                print_info(file=sys.stdout)
                progress_printer = RunningUI()
                progress_printer.max_threads = max_threads
                progress_printer.start_time = ui.start_time
                progress_printer.refresh(nodegraph)
                progress_printer.flush()
            elif character in "hH":
                if help_printed:
                    continue

                help_printed = True
                print_info("""
Commands:
  Key   Function
  h     Prints this message.
  l     Lists the currently runnning nodes.
  +     Increases the maximum number of threads by one.
  -     Decreases the maximum number of threads by one; already running tasks
        are NOT terminated if the number of threads currently used exceeds the
        resulting maximum.
""", file=sys.stdout)
            else:
                continue

        if max_threads != old_max_threads:
            print_debug("Maximum number of threads changed from %i to %i."
                        % (old_max_threads, max_threads), file=sys.stdout)

        return max_threads

    @classmethod
    def poll_stdin(cls):
        return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])


class BaseUI(object):
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
        self.states = []
        self.threads = 0
        self.max_threads = 0
        self.start_time = None
        self._end_time = None
        self._updated = True

    def flush(self):
        """Called by the user of the UI to ensure that the UI to print
        the current state of the pipeline / changes to pipeline / etc.

        Returns true if node-states have changed since last update.
        """
        if self._updated:
            self._updated = False
            return True
        return False

    def finalize(self):
        """Called by the pipeline at the termination of a run. By default,
        this function prints the location of the log-file if one was created
        during the run (e.g. if there were errors), and a summary of all nodes.
        """
        runtime = (self._end_time or 0) - (self.start_time or 0)

        if self.states[self.ERROR]:
            print_err("Done; but errors were detected ...")
        else:
            print_info("Done ...")

        print_info()
        rows = [("  Number of nodes:", sum(self.states)),
                ("  Number of done nodes:", self.states[self.DONE]),
                ("  Number of runable nodes:", self.states[self.RUNABLE]),
                ("  Number of queued nodes:", self.states[self.QUEUED]),
                ("  Number of outdated nodes:", self.states[self.OUTDATED]),
                ("  Number of failed nodes:", self.states[self.ERROR]),
                ("  Pipeline runtime:", _fmt_runtime(runtime))]

        for line in text.padded_table(rows):
            print_info(line)

        print_info("\nUse --list-output-files to view status of output files.")

        logfile = paleomix.logger.get_logfile()
        if logfile:
            print_debug("Log-file located at %r" % (logfile,))

        print_info()

    def refresh(self, nodegraph):
        """Called when the nodegraph has refreshed, causing state-counts
        to be recalculated."""
        self._updated = True
        self.states, self.threads \
            = self._count_states(nodegraph, nodegraph.iterflat())

    def state_changed(self, node, old_state, new_state, _is_primary):
        """Observer function for NodeGraph; counts states of nodes."""
        self._updated = True

        self.states[old_state] -= 1
        self.states[new_state] += 1
        if old_state == self.RUNNING:
            self.threads -= node.threads
        elif new_state == self.RUNNING:
            self.threads += node.threads

        if self.start_time is None:
            self.start_time = time.time()
        self._end_time = time.time()

    @classmethod
    def _count_states(cls, nodegraph, nodes):
        """Counts the number of each state observed for a set of nodes, and
        returns these as a list, as well as the estimated number of threads
        being used by running nodes."""
        states = [0] * nodegraph.NUMBER_OF_STATES
        threads = 0

        for node in nodes:
            state = nodegraph.get_node_state(node)

            states[state] += 1
            if state == nodegraph.RUNNING:
                threads += node.threads

        return states, threads

    def _describe_state(self):
        """ TODO """
        runtime = 0
        if self.start_time is not None:
            runtime = time.time() - self.start_time

        fields = [datetime.datetime.now().strftime("%T"),
                  ' Running ', str(self.states[self.RUNNING]), ' ',
                  ('task ' if self.states[self.RUNNING] == 1 else 'tasks '),
                  'using ~%i of max %i threads; ' % (self.threads,
                                                     self.max_threads)]

        if self.states[self.OUTDATED]:
            fields.append('%i outdated, ' % (self.states[self.OUTDATED],))

        if self.states[self.ERROR]:
            fields.append('%i failed, ' % (self.states[self.ERROR],))

        fields.extend(('%i done of %i tasks' % (self.states[self.DONE],
                                                sum(self.states),),
                       ' in ', _fmt_runtime(runtime),
                       '; press \'h\' for help.'))

        return ''.join(fields)

    DONE = paleomix.nodegraph.NodeGraph.DONE
    RUNNING = paleomix.nodegraph.NodeGraph.RUNNING
    RUNABLE = paleomix.nodegraph.NodeGraph.RUNABLE
    QUEUED = paleomix.nodegraph.NodeGraph.QUEUED
    OUTDATED = paleomix.nodegraph.NodeGraph.OUTDATED
    ERROR = paleomix.nodegraph.NodeGraph.ERROR


class RunningUI(BaseUI):
    """Prints a summary, and the list of running nodes every
    time 'flush' is called."""

    def __init__(self):
        BaseUI.__init__(self)
        self._running_nodes = []

    def flush(self):
        """See BaseUI.flush."""
        if BaseUI.flush(self) and self._running_nodes:
            self._print_header()
            for node in sorted(map(str, self._running_nodes)):
                print_info("  - %s" % (node,), file=sys.stdout)
            print_info(file=sys.stdout)

    def refresh(self, nodegraph):
        """See BaseUI.refresh."""
        BaseUI.refresh(self, nodegraph)
        self._running_nodes = []
        for node in nodegraph.iterflat():
            if nodegraph.get_node_state(node) == self.RUNNING:
                self._running_nodes.append(node)

    def state_changed(self, node, old_state, new_state, is_primary):
        """See BaseUI.state_changed."""
        BaseUI.state_changed(self, node, old_state, new_state, is_primary)

        if old_state == self.RUNNING:
            self._running_nodes.remove(node)
        elif new_state == self.RUNNING:
            self._running_nodes.append(node)

    def _print_header(self):
        print_msg('\n%s' % (self._describe_state(),), file=sys.stdout)

        logfile = paleomix.logger.get_logfile()
        if logfile:
            print_debug("  Log-file located at %r" % (logfile,),
                        file=sys.stdout)


class ProgressUI(BaseUI):
    """Progress based UI: Prints nodes when they start running; they finish
    running; or when they fail running. Changes to state resulting from the
    above is not printed. Every 20th update is followed by a summary of the
    current total progress when flush is called."""

    # Print a summery of the current state very N events
    _SUMMARY_EVERY = 20

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
            self._print_state(node, new_state)

            self._refresh_count -= 1
            if new_state == self.ERROR:
                self._refresh_count = 0

    def flush(self):
        """See BaseUI.flush."""
        if (self._refresh_count <= 0):
            self._refresh_count = ProgressUI._SUMMARY_EVERY
            self._print_summary()

    def _print_summary(self):
        """Prints a summary of the pipeline progress."""
        print_msg()
        print_msg(self._describe_state(), file=sys.stdout)

        logfile = paleomix.logger.get_logfile()
        if logfile:
            print_debug("Log-file located at %r" % (logfile,), file=sys.stdout)

    def _print_state(self, node, new_state):
        state_label, print_func = self._DESCRIPTIONS[new_state]
        if new_state == self.RUNNING:
            self._runtimes[node] = time.time()
        elif new_state in (self.RUNNING, self.DONE, self.ERROR):
            state_label = "%s (%s)" % (state_label, self._get_runtime(node))

        time_label = datetime.datetime.now().strftime("%T")
        print_func("%s %s: %s" % (time_label, state_label, node),
                   file=sys.stdout)

    def _get_runtime(self, node):
        current_time = time.time()
        runtime = current_time - self._runtimes.pop(node, current_time)
        return _fmt_runtime(runtime)

    _DESCRIPTIONS = {
        BaseUI.DONE: ("Finished", print_disabled),
        BaseUI.RUNNING: ("Started", print_info),
        BaseUI.ERROR: ("Failed", print_err),
    }


class SummaryUI(BaseUI):
    def __init__(self):
        self._max_len = 0
        self._new_error = False
        BaseUI.__init__(self)

    def state_changed(self, node, old_state, new_state, is_primary):
        BaseUI.state_changed(self, node, old_state, new_state, is_primary)
        self._new_error |= (new_state == self.ERROR and is_primary)

    def flush(self):
        if BaseUI.flush(self):
            description = self._describe_state()

            self._max_len = max(len(description), self._max_len)
            print_msg("\r%s" % (description.ljust(self._max_len),), end="",
                      file=sys.stdout)

            logfile = paleomix.logger.get_logfile()
            if logfile and self._new_error:
                print_debug("\nLog-file located at %r" % (logfile,),
                            file=sys.stdout)
                self._new_error = False
            sys.stdout.flush()

    def finalize(self):
        print_msg(file=sys.stdout)
        BaseUI.finalize(self)


def _fmt_runtime(runtime):
    runtime = int(round(runtime))

    if runtime >= 3600:
        fmt = "{hours}:{mins:02}:{secs:02}s"
    elif runtime >= 60:
        fmt = "{mins}:{secs:02}s"
    else:
        fmt = "{secs}s"

    return fmt.format(hours=int(runtime) // 3600,
                      mins=(int(runtime) // 60) % 60,
                      secs=(runtime % 60))


# No longer provided
VerboseUI = RunningUI
QuietUI = RunningUI

# Different types of UIs
UI_TYPES = {
    "Verbose": VerboseUI,
    "Quiet": RunningUI,
    "Running": RunningUI,
    "Progress": ProgressUI,
    "Summary": SummaryUI,
}
