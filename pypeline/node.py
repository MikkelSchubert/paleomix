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
import os
import traceback

import pypeline.common.fileutils as fileutils
from pypeline.common.utilities import safe_coerce_to_tuple

# Imported from, in order to allow monkeypatching in tests
# FIXME: Make create_temp_dir a 'with' object
from os import rmdir
from pypeline.common.fileutils import create_temp_dir



class NodeError(RuntimeError):
    pass

class NodeUnhandledException(NodeError):
    """This exception is thrown by Node.run() if a non-NodeError exception 
    is raised in a subfunction (e.g. _setup, _run, or _teardown). The text
    for this exception will include both the original error message and a 
    stacktrace for that error."""
    pass



class Node(object):
    def __init__(self, description = None, threads = 1,
                 input_files = (), output_files = (),
                 executables = (), auxiliary_files = (),
                 requirements = (), subnodes = (), dependencies = ()):

        self.__description   = description
        self.input_files     = safe_coerce_to_tuple(input_files)
        self.output_files    = safe_coerce_to_tuple(output_files)
        self.executables     = safe_coerce_to_tuple(executables)
        self.auxiliary_files = safe_coerce_to_tuple(auxiliary_files)
        self.requirements    = safe_coerce_to_tuple(requirements)

        self.subnodes       = self._collect_nodes(subnodes, "Subnode")
        self.dependencies   = self._collect_nodes(dependencies, "Dependency")

        self.threads        = int(threads)


    @property
    def is_done(self):
        """Returns true if all subnodes of this node are done, and if all output
        files of this node exists (empty files are considered as valid). If the
        node doesn't produce output files, it is always considered done by. To
        change this behavior, override the 'is_done' property"""

        if not all(node.is_done for node in self.subnodes):
            return False
        elif fileutils.missing_files(self.output_files):
            return False

        return True


    @property
    def is_outdated(self):
        """Returns true if the output exists (is_done == True), but one or more
        of the intput files appear to have been changed since the creation of the
        output files (based on the timestamps). A node that lacks either input or
        output files is never considered outdated."""

        if not self.is_done:
            return False
        elif not (self.input_files and self.output_files):
            return False

        return fileutils.modified_after(self.input_files, self.output_files)


    def run(self, config):
        """Runs the node, by calling _setup, _run, and _teardown in that order.
        Prior to calling these functions, a temporary dir is created using the
        'temp_root' prefix from the config object. Both the config object and
        the temporary dir are passed to the above functions. The temporary
        dir is removed after _teardown is called, and all expected files 
        should have been removed/renamed at that point.

        Any non-NodeError exception raised in this function is wrapped in a
        NodeUnhandledException, which includes a full backtrace. This is needed
        to allow showing these in the main process."""
        
        try:
            temp = create_temp_dir(config.temp_root)
            
            self._setup(config, temp)
            self._run(config, temp)
            self._teardown(config, temp)

            rmdir(temp)
        except NodeError:
            raise
        except Exception:
            raise NodeUnhandledException(traceback.format_exc())


    def _setup(self, _config, _temp):
        """Is called prior to '_run()' by 'run()'. Any code used to copy/link files, 
        or other steps needed to ready the node for running may be carried out in this
        function. Checks that required input files exist, and raises an NodeError if 
        this is not the case."""
        if fileutils.missing_executables(self.executables):
            raise NodeError("Executable(s) does not exist for command: %s" \
                                % (self._command,))
        self._check_for_missing_files(self.input_files, "input")
        self._check_for_missing_files(self.auxiliary_files, "auxiliary")


    def _run(self, _config, _temp):
        pass

    def _teardown(self, _config, _temp):
        self._check_for_missing_files(self.output_files, "output")

    def __str__(self):
        """Returns the description passed to the constructor, or a default
        description if no description was passed to the constructor."""
        if self.__description:
            return self.__description
        return "<%s>" % (self.__class__.__name__,)


    def _collect_nodes(self, nodes, description):
        if nodes is None:
            return frozenset()

        nodes = frozenset(safe_coerce_to_tuple(nodes))
        bad_nodes = [node for node in nodes if not isinstance(node, Node)]

        if bad_nodes:
            bad_nodes = [repr(node) for node in bad_nodes]
            message = "%s-list contain non-Node objects:\n\t- Command: %s\n\t- Objects: %s" \
                % (description, self, "\n\t           ".join(bad_nodes))
            raise TypeError(message)

        return nodes


    def _check_for_missing_files(self, filenames, description):
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            message = "Missing %s files for command:\n\t- Command: %s\n\t- Files: %s" \
                % (description, self, "\n\t         ".join(missing_files))
            raise NodeError(message)


    @classmethod
    def _desc_files(cls, files):
        if len(files) == 1:
            return "'%s'" % tuple(files)
        else:
            paths = set(os.path.dirname(filename) for filename in files)
            if len(paths) == 1:
                return "%i file(s) in '%s'" % (len(files), paths.pop())
            else:
                return "%i file(s)" % (len(files),)




class CommandNode(Node):
    def __init__(self, command, description = None, threads = 1,
                 subnodes = (), dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      input_files  = command.input_files,
                      output_files = command.output_files,
                      auxiliary_files = command.auxiliary_files,
                      executables  = command.executables,
                      requirements = command.requirements,
                      threads      = threads,
                      subnodes     = subnodes,
                      dependencies = dependencies)

        self._command = command

    def _run(self, _config, temp):
        """Runs the command object provided in the constructor, and waits for it to 
        terminate. If any errors during the running of the command, this function
        raises a NodeError detailing the returned error-codes."""
        self._command.run(temp)

        return_codes = self._command.join()
        if any(return_codes):
            raise NodeError("Error(s) running '%s':\n\tReturn-codes: %s\n\tTemporary directory: '%s'" \
                             % (self._command, return_codes, temp))


    def _teardown(self, config, temp):
        self._command.commit(temp)

        Node._teardown(self, config, temp)



class MetaNode(Node):
    """A MetaNode is a simplified node, which only serves the purpose of aggregating
    a set of subnodes. It does not carry out any task itself (run() does nothing),
    and is marked as done when all its subnodes / dependencies are completed."""

    def __init__(self, description = None, subnodes = (), dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      subnodes     = subnodes,
                      dependencies = dependencies)

    @property
    def is_done(self):
        raise RuntimeError("Called 'is_done' on MetaNode")

    @property
    def is_outdated(self):
        raise RuntimeError("Called 'is_outdated' on MetaNode")

    def run(self, config):
        raise RuntimeError("Called 'run' on MetaNode")
