import os
import subprocess


import ui
import fileutils


SINGLE_THREAD = 1


class Node:
    def __init__(self, description = None, command = [], dependencies = []):
        self.__description = description
        self.__command = command
        self.dependencies = dependencies


    def get_commands(self):
        """Returns a list of sets of commands (AtomicCmd or AtomicSet),
        each of which is run to completion, before the next set is requested. 
        By default, only the commands set in the constructor are returned. 
        Functions  doing inter-command processing may overload this function, 
        rather than passing values to the constructor.

        If this function is overloaded to add inter-command processing,
        the 'out_exists' function should also be overloaded."""
        yield self.__command


    def output_exists(self):
        for command in self.get_commands():
            if command.missing_output_files():
                return False
        return True


    def __str__(self):
        if self.__description:
            return self.__description
        return "<%s>" % (self.__class__.__name__,)


    def run(self):
        ui.print_info("%s: Running ..." % ((self,)))

        temp = fileutils.create_temp_dir()
        for command in self.get_commands():
            if not self._check_prerequisites(command):
                os.rmdir(temp)
                return False
            
            command.run(temp)
            
            if not self._wait_and_commit(command):
                ui.print_warn("%s: Failed ..." % (self,))
                return False

        os.rmdir(temp)        
        ui.print_msg("%s: Finished" % ((self,)))
        return True


    def _check_prerequisites(self, command):
        missing_input = command.missing_input_files()
        if missing_input:
            ui.print_warn("%s: Missing input files for command: %s -> %s" \
                              % (self, command, missing_input))
            return False

        return True

    def _wait_and_commit(self, command):
        return_codes = command.wait()
        if any(return_codes):
            ui.print_err("%s: Error(s) running '%s': Return-codes %s" \
                             % (self, command, return_codes))
            return False
        
        missing_files = command.missing_temp_files()
        if missing_files:
            ui.print_err("%s: Expected temp files are missings '%s'" \
                             % (self, command, missing_files))
            return False
    
        command.commit()
        return True
