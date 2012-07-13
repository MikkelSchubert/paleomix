import os
import traceback

import fileutils

from pylib.utilities import safe_coerce_to_tuple




class NodeError(RuntimeError):
    def __init__(self, *vargs):
        RuntimeError.__init__(self, *vargs)


class NodeUnhandledException(NodeError):
    def __init__(self, *vargs):
        NodeError.__init__(self, *vargs)



class Node(object):
    def __init__(self, description = None, 
                 input_files = (), output_files = (), 
                 subnodes = (), dependencies = ()):

        self.__description  = description
        self.__input_files  = safe_coerce_to_tuple(input_files)
        self.__output_files = safe_coerce_to_tuple(output_files)

        self.subnodes       = self._collect_nodes(subnodes, "Subnode")
        self.dependencies   = self._collect_nodes(dependencies, "Dependency")


    @property
    def is_done(self):
        for node in self.subnodes:
            if not node.is_done:
                return False

        if fileutils.missing_files(self.__output_files):
            return False

        return True


    @property
    def is_outdated(self):
        if not self.is_done:
            return False
        elif not (self.__input_files and self.__output_files):
            return False

        return fileutils.modified_after(self.__input_files, self.__output_files)


    @property
    def input_files(self):
        return self.__input_files


    @property
    def output_files(self):
        return self.__output_files


    def run(self, config):
        try:
            temp = fileutils.create_temp_dir(config.temp_root)
            
            self._setup(config, temp)
            self._run(config, temp)
            self._teardown(config, temp)

            os.rmdir(temp)
        except NodeError:
            raise
        except Exception:
            raise NodeUnhandledException(traceback.format_exc())


    def _setup(self, _config, _temp):
        self._check_for_missing_files(self.__input_files, "input")

    def _run(self, _config, _temp):
        pass

    def _teardown(self, _config, _temp):
        self._check_for_missing_files(self.__output_files, "output")

    def __str__(self):
        if self.__description:
            return self.__description
        return "<%s>" % (self.__class__.__name__,)


    def _collect_nodes(self, nodes, description):
        nodes = frozenset(safe_coerce_to_tuple(nodes))
        bad_nodes = [node for node in nodes if not isinstance(node, Node)]

        if bad_nodes:
            message = "%s-list contain non-Node objects:\n\t- Command: %s\n\t- Objects: %s" \
                % (description, self, "\n\t           ".join(map(repr, bad_nodes)))
            raise TypeError(message)

        return nodes


    def _check_for_missing_files(self, filenames, description):
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            message = "Missing %s files for command:\n\t- Command: %s\n\t- Files: %s" \
                % (description, self, "\n\t         ".join(missing_files))
            raise NodeError(message)





class CommandNode(Node):
    def __init__(self, command, description = None, 
                 subnodes = (), dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      input_files  = command.input_files(),
                      output_files = command.output_files(),
                      subnodes     = subnodes,
                      dependencies = dependencies)

        self._command = command
            

    def _setup(self, config, temp):
        if fileutils.missing_executables(self._command.executables()):
            raise NodeError("Executable(s) does not exist for command: %s" \
                                % (self._command, ))

        Node._setup(self, config, temp)


    def _run(self, _config, temp):
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
    and is marked as done when all its subnodes are completed."""

    def __init__(self, description = None, subnodes = (), dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      subnodes     = subnodes,
                      dependencies = dependencies)

    @property
    def is_done(self):
        for node in self.subnodes:
            if not node.is_done:
                return False

        return True


    def run(self, config):
        pass
