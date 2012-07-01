import os
import traceback

import ui
import fileutils


SINGLE_THREAD = 1

def _validate_nodes(nodes, description):
    nodes = tuple(nodes)
    for dependency in nodes:
        if not isinstance(dependency, Node):
            raise NodeError("%s is not a Node: %s" \
                                % (description, dependency))
    return nodes


class NodeError(RuntimeError):
    def __init__(self, *vargs):
      RuntimeError.__init__(self, *vargs)


class Node:
    def __init__(self, description = None, dependencies = ()):
        self.__description = description
        self.dependencies = _validate_nodes(dependencies, "Dependency")
    
    def output_exists(self):
        assert False

    def run(self, config):
        try:
            temp = fileutils.create_temp_dir(config.temp_root)
            
            if not self._setup(config, temp):
                return False
            elif not self._run(config, temp):
                return False
            elif not self._teardown(config, temp):
                return False

            os.rmdir(temp)
            
            return True
        except Exception:
            raise NodeError(traceback.format_exc())


    def _setup(self, _config, _temp):
        return True

    def _run(self, _config, _temp):
        return True

    def _teardown(self, _config, _temp):
        return True


    def __str__(self):
        if self.__description:
            return self.__description
        return "<%s>" % (self.__class__.__name__,)




class SimpleNode(Node):
    def __init__(self, command = None, description = None, optional_files = (), dependencies = ()):
        Node.__init__(self, description, dependencies)

        self._command = command
        self._optional_files = frozenset(optional_files)

    
    def output_exists(self):
        assert self._command
        missing_files = self._command.missing_output_files()
        for optional in self._optional_files:
            if optional in missing_files:
                missing_files.remove(optional)

        return not missing_files


    def _setup(self, _config, _temp):
        if not self._command.executable_exists():
            ui.print_err("%s: Executable does not exist for command: %s" \
                              % (self, self._command))
            return False
            
        missing_input = self._command.missing_input_files()
        if missing_input:
            ui.print_warn("%s: Missing input files for command: %s -> %s" \
                              % (self, self._command, missing_input))
            return False

        return True


    def _run(self, _config, temp):
        self._command.run(temp)

        return_codes = self._command.wait()
        if any(return_codes):
            ui.print_err("%s: Error(s) running '%s': Return-codes %s" \
                             % (self, self._command, return_codes))
            return False

        return True


    def _teardown(self, _config, _temp):
        missing_files = self._command.missing_temp_files()
        if missing_files:
            ui.print_err("%s: Expected temp files are missings for command '%s': %s" \
                             % (self, self._command, missing_files))
            return False
    
        return self._command.commit()
