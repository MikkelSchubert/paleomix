import os
import traceback

import fileutils



class NodeError(RuntimeError):
    def __init__(self, *vargs):
        RuntimeError.__init__(self, *vargs)


class NodeUnhandledException(RuntimeError):
    def __init__(self, *vargs):
        RuntimeError.__init__(self, *vargs)



class Node(object):
    def __init__(self, description = None, input_files = (), output_files = (), dependencies = ()):
        self.__description  = description
        self.__input_files  = tuple(input_files)
        self.__output_files = tuple(output_files)

        if isinstance(dependencies, Node):
            dependencies = (dependencies,)

        try:
            self.subnodes = tuple(dependencies)
        except TypeError:
            raise NodeError("Non Node/sequence passed as dependencies for '%s': %s" \
                    % (self, repr(dependencies)))

        for subnode in self.subnodes:
            if not isinstance(subnode, Node):
                raise NodeError("Subnode is not a Node: %s" % subnode)


    @property
    def is_done(self):
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


    def run(self, config):
        try:
            temp = fileutils.create_temp_dir(config.temp_root)
            
            self._setup(config, temp)
            self._run(config, temp)
            self._teardown(config, temp)

            os.rmdir(temp)
            
            return True
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


    def _check_for_missing_files(self, filenames, description):
        missing_files = fileutils.missing_files(filenames)
        if missing_files:
            message = "Missing %s files for command:\n\t- Command: %s\n\t- Files: %s" \
                % (description, self, "\n\t         ".join(missing_files))
            raise NodeError(message)





class CommandNode(Node):
    def __init__(self, command, description = None, dependencies = ()):
        Node.__init__(self, 
                      description  = description,
                      input_files  = command.input_files(),
                      output_files = command.output_files(),
                      dependencies = dependencies)

        self._command = command
            

    def _setup(self, config, temp):
        if fileutils.missing_executables(self._command.executables()):
            raise NodeError("Executable(s) does not exist for command: %s" \
                                % (self._command, ))

        Node._setup(self, config, temp)


    def _run(self, _config, temp):
        self._command.run(temp)

        return_codes = self._command.wait()
        if any(return_codes):
            raise NodeError("Error(s) running '%s':\n\tReturn-codes: %s\n\tTemporary directory: '%s'" \
                             % (self._command, return_codes, temp))


        def _teardown(self, config, temp):
        self._check_for_missing_files(self._command.temp_files(), "")
        self._command.commit()

        Node._teardown(self, config, temp)



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


    def run(self, config):
        pass
