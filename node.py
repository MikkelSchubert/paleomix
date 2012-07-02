import os
import traceback

import ui
import fileutils



class NodeError(RuntimeError):
    def __init__(self, *vargs):
        RuntimeError.__init__(self, *vargs)


class Node(object):
    def __init__(self, description = None, input_files = (), output_files = (), dependencies = ()):
        self.__description  = description
        self.__input_files  = tuple(input_files)
        self.__output_files = tuple(output_files)

        self.subnodes = tuple(dependencies)
        for subnode in self.subnodes:
            if not isinstance(subnode, Node):
                raise NodeError("Subnode is not a Node: %s" % subnode)


    @property
    def is_done(self):
        if fileutils.missing_files(self.__output_files):
            return False
        elif self.is_outdated:
            return False

        return True


    @property
    def is_runable(self):
        for subnode in self.subnodes:
            if not subnode.is_done:
                return False

        return True


    @property
    def is_outdated(self):
        if fileutils.missing_files(self.__output_files):
            return False

        for node in self.subnodes:
            if node.is_outdated or not node.is_done:
                return True

        if not (self.__input_files and self.__output_files):
            return False

        return fileutils.modified_after(self.__input_files, self.__output_files)


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




class CommandNode(Node):
    def __init__(self, command, description = None, dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      input_files  = command.input_files(),
                      output_files = command.output_files(),
                      dependencies = dependencies)

        self._command = command
            

    def _setup(self, _config, _temp):
        if not fileutils.missing_files(self._command.executables()):
            ui.print_err("%s: Executable(s) does not exist for command: %s" \
                              % (self, self._command))
            return False
            
        missing_input = fileutils.missing_files(self._command.input_files())
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
            ui.print_err("\t - Ran in directory '%s'." % temp)
            return False

        return True


    def _teardown(self, _config, _temp):
        missing_files = fileutils.missing_files(self._command.temp_files())
        if missing_files:
            ui.print_err("%s: Expected temp files are missings for command '%s': %s" \
                             % (self, self._command, missing_files))
            return False
    
        return self._command.commit()



class MetaNode(Node):
    def __init__(self, description = None, subnodes = ()):
        Node.__init__(self, 
                      description  = description,
                      dependencies = tuple(subnodes))

    @property
    def is_done(self):
        for node in self.subnodes:
            if not node.is_done:
                return False

        return True


    def run(self):
        return True
