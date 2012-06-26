import os

import ui
import fileutils


SINGLE_THREAD = 1


class Node:
    def __init__(self, description = None, command = None, optional_files = (), dependencies = ()):
        self.__description = description
        self.__command = command
        self.__optional_files = optional_files
        self.dependencies = tuple(dependencies)

    def output_exists(self):
        assert self.__command
        missing_files = self.__command.missing_output_files()
        for optional in self.__optional_files:
            if optional in missing_files:
                missing_files.remove(optional)

        return not missing_files

    def run(self, config):
        temp = fileutils.create_temp_dir(config.temp_root)
        if not self._check_prerequisites(self.__command):
            os.rmdir(temp)
            return False
            
        self.__command.run(temp)
            
        if not self._wait_for_command(self.__command):
            return False
        elif not self._commit_command(self.__command):
            return False

        os.rmdir(temp)
        return True


    def __str__(self):
        if self.__description:
            return self.__description
        return "<%s>" % (self.__class__.__name__,)


    def _check_prerequisites(self, command):
        if not command.executable_exists():
            ui.print_err("%s: Executable does not exist for command: %s" \
                              % (self, command))
            return False
            
        missing_input = command.missing_input_files()
        if missing_input:
            ui.print_warn("%s: Missing input files for command: %s -> %s" \
                              % (self, command, missing_input))
            return False

        return True


    def _wait_for_command(self, command):
        return_codes = command.wait()
        if any(return_codes):
            ui.print_err("%s: Error(s) running '%s': Return-codes %s" \
                             % (self, command, return_codes))
            return False

        return True


    def _commit_command(self, command):
        missing_files = command.missing_temp_files()
        if missing_files:
            ui.print_err("%s: Expected temp files are missings for command '%s': %s" \
                             % (self, command, missing_files))
            return False
    
        command.commit()
        return True
