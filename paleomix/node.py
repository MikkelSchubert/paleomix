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
import logging
import os
import shutil
import sys
import traceback
from pathlib import Path

import paleomix
import paleomix.common.fileutils as fileutils
from paleomix.atomiccmd.command import CmdError
from paleomix.common.utilities import safe_coerce_to_frozenset


class NodeError(RuntimeError):
    pass


class CmdNodeError(NodeError):
    pass


class NodeUnhandledException(NodeError):
    """This exception is thrown by Node.run() if a non-NodeError exception
    is raised in a subfunction (e.g. _setup, _run, or _teardown). The text
    for this exception will include both the original error message and a
    stacktrace for that error."""

    pass


class Node:
    def __init__(
        self,
        description=None,
        threads=1,
        input_files=(),
        output_files=(),
        executables=(),
        auxiliary_files=(),
        requirements=(),
        dependencies=(),
    ):
        if not isinstance(description, _DESC_TYPES):
            raise TypeError(
                "'description' must be None or a string, not %r"
                % (description.__class__.__name__,)
            )

        self.__description = description
        self.input_files = self._validate_files(input_files)
        self.output_files = self._validate_files(output_files)
        self.executables = self._validate_files(executables)
        self.auxiliary_files = self._validate_files(auxiliary_files)
        self.requirements = self._validate_requirements(requirements)

        self.threads = self._validate_nthreads(threads)
        self.dependencies = self._collect_nodes(dependencies)

        # If there are no input files, the node cannot be re-run based on
        # changes to the input, and nodes with output but no input are not
        # expected based on current usage.
        if not self.input_files and self.output_files:
            raise NodeError("Node not dependent upon input files: %s" % self)

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

        temp = None

        try:
            # Generate directory name and create dir at temp_root
            temp = self._create_temp_dir(config)

            self._setup(config, temp)
            self._run(config, temp)
            self._teardown(config, temp)
            self._remove_temp_dir(temp)
        except NodeError as error:
            self._write_error_log(temp, error)
            raise NodeError(
                "Error(s) running Node:\n\tTemporary directory: %s\n\n%s"
                % (repr(temp), error)
            )

        except Exception as error:
            self._write_error_log(temp, error)
            raise NodeUnhandledException(
                "Error(s) running Node:\n\tTemporary directory: %s\n\n%s"
                % (repr(temp), traceback.format_exc())
            )

    def _create_temp_dir(self, config):
        """Called by 'run' in order to create a temporary folder.

        'config' is expected to have a property .temp_root under
        which the temporary folder is created."""
        return fileutils.create_temp_dir(config.temp_root)

    def _remove_temp_dir(self, temp):
        """Called by 'run' in order to remove an (now) empty temporary folder."""
        temp = fileutils.fspath(temp)
        log = logging.getLogger(__name__)
        for filename in self._collect_files(temp):
            log.warning(
                "Unexpected file in temporary directory: %r",
                os.path.join(temp, filename),
            )

        shutil.rmtree(temp)

    def _setup(self, _config, _temp):
        """Is called prior to '_run()' by 'run()'. Any code used to copy/link files,
        or other steps needed to ready the node for running may be carried out in this
        function. Checks that required input files exist, and raises an NodeError if
        this is not the case."""
        if fileutils.missing_executables(self.executables):
            raise NodeError("Executable(s) does not exist for node: %s" % (self,))

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
        return repr(self)

    def __getstate__(self):
        """Called by pickle/cPickle to determine what to pickle; this is
        overridden to avoid pickling of requirements, dependencies, which would
        otherwise greatly inflate the amount of information that needs to be
        pickled."""
        obj_dict = self.__dict__.copy()
        obj_dict["requirements"] = ()
        obj_dict["dependencies"] = ()
        return obj_dict

    def _write_error_log(self, temp, error):
        if not (temp and os.path.isdir(temp)):
            return

        def _fmt(values):
            return "\n                   ".join(sorted(values))

        message = [
            "PALEOMIX         = v%s" % (paleomix.__version__,),
            "Command          = %r" % (" ".join(sys.argv),),
            "CWD              = %r" % (os.getcwd(),),
            "PATH             = %r" % (os.environ.get("PATH", ""),),
            "Node             = %s" % (str(self),),
            "Threads          = %i" % (self.threads,),
            "Input files      = %s" % (_fmt(self.input_files),),
            "Output files     = %s" % (_fmt(self.output_files),),
            "Auxiliary files  = %s" % (_fmt(self.auxiliary_files),),
            "Executables      = %s" % (_fmt(self.executables),),
            "",
            "Errors =\n%s\n" % (error,),
        ]
        message = "\n".join(message)

        try:
            with open(os.path.join(temp, "pipe.errors"), "w") as handle:
                handle.write(message)
        except OSError as oserror:
            sys.stderr.write("ERROR: Could not write failure log: %s\n" % (oserror,))

    def _collect_nodes(self, nodes):
        if nodes is None:
            return frozenset()

        nodes = safe_coerce_to_frozenset(nodes)
        bad_nodes = [node for node in nodes if not isinstance(node, Node)]

        if bad_nodes:
            bad_nodes = [repr(node) for node in bad_nodes]
            message = (
                "Dependency-list contain non-Node objects:\n"
                "\t- Command: %s\n\t- Objects: %s"
                % (self, "\n\t           ".join(bad_nodes))
            )
            raise TypeError(message)

        return nodes

    def _check_for_missing_files(self, filenames, description):
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            message = (
                "Missing %s files for command:\n\t- Command: %s\n\t- Files: %s"
                % (description, self, "\n\t         ".join(missing_files))
            )
            raise NodeError(message)

    @classmethod
    def _validate_requirements(cls, requirements):
        requirements = safe_coerce_to_frozenset(requirements)
        for requirement in requirements:
            if not callable(requirement):
                raise TypeError(
                    "'requirements' must be callable, not %r" % (type(requirement),)
                )
        return requirements

    @classmethod
    def _validate_files(cls, files):
        return frozenset(fileutils.validate_filenames(files))

    @classmethod
    def _validate_nthreads(cls, threads):
        if not isinstance(threads, int):
            raise TypeError(
                "'threads' must be a positive integer, not a %s" % (type(threads),)
            )
        elif threads < 1:
            raise ValueError(
                "'threads' must be a positive integer, not %i" % (threads,)
            )
        return threads

    @staticmethod
    def _collect_files(root):
        root = fileutils.fspath(root)

        def _walk_dir(path):
            for entry in os.scandir(path):
                if entry.is_file():
                    yield str(Path(entry.path).relative_to(root))
                elif entry.is_dir():
                    yield from _walk_dir(entry.path)

        yield from _walk_dir(root)


class CommandNode(Node):
    def __init__(self, command, description=None, threads=1, dependencies=()):
        Node.__init__(
            self,
            description=description,
            input_files=command.input_files,
            output_files=command.output_files,
            auxiliary_files=command.auxiliary_files,
            executables=command.executables,
            requirements=command.requirements,
            threads=threads,
            dependencies=dependencies,
        )

        self._command = command

    def _run(self, _config, temp):
        """Runs the command object provided in the constructor, and waits for it to
        terminate. If any errors during the running of the command, this function
        raises a NodeError detailing the returned error-codes."""
        try:
            self._command.run(temp)
        except CmdError as error:
            raise CmdNodeError("%s\n\n%s" % (str(self._command), error))

        return_codes = self._command.join()
        if any(return_codes):
            raise CmdNodeError(str(self._command))

    def _teardown(self, config, temp):
        required_files = self._command.expected_temp_files
        current_files = set(self._collect_files(temp))

        missing_files = required_files - current_files
        if missing_files:
            raise CmdNodeError(
                (
                    "Error running Node, required files not created:\n"
                    "Temporary directory: %r\n"
                    "\tRequired files missing from temporary directory:\n\t    - %s"
                )
                % (temp, "\n\t    - ".join(sorted(map(repr, missing_files))))
            )

        self._command.commit(temp)

        Node._teardown(self, config, temp)


# Types that are allowed for the 'description' property
_DESC_TYPES = (str, type(None))
