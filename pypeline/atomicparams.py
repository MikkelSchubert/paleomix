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
"""Tools for passing CLI options to AtomicCmds used by Nodes.

The module contains 1 class and 2 decorators, which may be used in conjunction
to create Node classes for which the call carried out by AtomicCmds may be 
modified by the end user, without explicit support added to the init function
of the class. The basic outline of such a class is as follows:


class ExampleNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(self, ...):
        # Use passed parameters to create AtomicParams obj
        params = AtomicParams(...)
        params.set_parameter(...)
        
        # Return dictionary of AtomicParams objects and any
        # additional parameters required to run the Node. 
        return {"command"   : params,
                "example" : ...}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        # Create AtomicCmd object using (potentially tweaked) parameters
        command = parameters.command.finalize()

        # Do something with a parameter passed to customize
        description = "<ExampleNode: %s>" % parameters.example
        
        CommandNode.__init__(command     = command,
                             description = description,
                             ...)

This class can then be used in two ways:
1) Without doing any explicit modifications to the CLI calls:
>> node = ExampleNode(...)

2) Retrieving and tweaking AtomicParams before creating the Node:
>> params = ExampleNode.customize(...)
>> params.command.set_option(...)
>> node = params.build_node()

"""
import os
import random
import inspect
import itertools
import collections

from pypeline.atomiccmd import AtomicCmd
from pypeline.common.utilities import safe_coerce_to_tuple



class ParameterError(RuntimeError):
    pass



class AtomicParams:
    """AtomicParams is a class used to construct a system call by a Node for use
    with its AtomicCmd(s). This allows the user of a Node to modify the behavior
    of the called programs using some CLI parameters, without explicit support 
    for these in the Node interface. Some limitations are in place, to help 
    catch cases where overwriting or adding a flag would break the Node.

    Parameters are divided into two classes; singleton and non-singleton.
     - Singleton parameters are those set exactly once on the command-line, 
       with subsequent calls overwriting the previous value. The function
       'set_parameter' is used to set such parameters.
     - Non-singleton parameters are those which may be set any number of 
       times, with all values specified being used by the program. Such
       parameters are added using the 'push_parameter' function.

    In addition to parameters, path may be assigned to keywords (for use as
    values for parameters), corresponding to their use in AtomicCmd.
    """

    def __init__(self, call):
        """Initialize AtomicParams object.

        The 'call' parameter specifies the basic call, and typically contains
        just the name of the executable. It may also be used to invoke a script
        or similar using an executable, while allowing the user to pass 
        parameters to the script (see e.g. AtomicJavaParams)."""
        
        self._call    = safe_coerce_to_tuple(call)
        self._fields  = []
        self._params  = []
        self._paths   = {}
        self._cmd_object = None


    def push_parameter(self, key, value = None, sep = None, fixed = True):
        """Appends a parameter that may be specified any number of times.

        This function is meant for parameters that may be specified any number of times,
        as exemplified by the "I" parameter used by many PicardTools (i.e. a non-singleton
        parameter), and a ParamterError exception will be raised if the parameter has 
        previously been set using 'set_parameter' (i.e. a singleton parameter).

        If no value has been specified (is None), only the key will added to the system-
        call. If a value and a sep(erator) has been specified, the key and the value will
        be joined into a single string before being added to the system-call. If fixed is 
        true, this particular parameter cannot be removed using 'pop_parameter'.
        
        Example:
        >> params = ...
        >> params.push_paramter("I", "%(IN_FILE_1)s", sep = "=")
        >> params.push_paramter("I", "%(IN_FILE_2)s", sep = "=")
        >> params.call
        [..., "I=%(IN_FILE_1)s", "I=%(IN_FILE_2)s"]

        Use keywords (as above) and 'set_paths' to ensure that input/output files are tracked!
        """
        if self._cmd_object:
            raise ParameterError("Parameters have already been finalized")
        elif any(opt["Singleton"] for opt in self._params if (opt["Key"] == key)):
            raise ParameterError("Attempted to push parameters previously marked as singleton: %s" % key)

        self._params.append({"Key"    : key,
                             "Value"  : value,
                             "Sep"    : sep,
                             "Fixed"  : fixed,
                             "Singleton" : False})


    def set_parameter(self, key, value = None, sep = None, fixed = True):
        """Sets or overwrites a parameter that may be specified at most once.


        This function is meant for parameters that should only be specified once, as
        is the typical expectation of CLI swithces (e.g. RAxML '-f'). A ParameterError
        will be raised if the option has previously been added using 'push_parameter'.

        If no value has been specified (is None), only the key will added to the system-
        call. If a value and a sep(erator) has been specified, the key and the value will
        be joined into a single string before being added to the system-call. If fixed is 
        true, this particular parameter cannot be removed using 'pop_parameter', or 
        overwritten by subsequent calls to 'set_parameter'."""
        if self._cmd_object:
            raise ParameterError("Parameters have already been finalized")
        elif any(not opt["Singleton"] for opt in self._params if (opt["Key"] == key)):
            raise ParameterError("Attempted to overwrite non-singleton parameter: %s" % key)

        param = {"Key" : key, "Value"  : value, "Sep" : sep, "Fixed"  : fixed, "Singleton" : True}
        for o_param in self._params:
            if o_param["Key"] == key:
                if o_param["Fixed"]:
                    raise ParameterError("Attempted to overwrite fixed singleton parameter: %s" % key)
                o_param.update(param)
                return

        self._params.append(param)


    def pop_parameter(self, key):
        if self._cmd_object:
            raise ParameterError("Parameters have already been finalized")

        for parameter in reversed(self._params):
            if parameter["Key"] == key:
                if parameter["Fixed"]:
                    raise ParameterError("Attempted to pop fixed parameter: %s" % key)
                self._params.remove(parameter)
                return

        raise KeyError("Could not pop parameter: %s" % key)


    def set_paths(self, key = None, value = None, **kwargs):
        if self._cmd_object:
            raise ParameterError("Parameters have already been finalized")
        elif key and value:
            kwargs[key] = value
        elif key or value:
            raise ParameterError("Either neither or both 'key' and 'value' parameters must be specified.")
        
        for (key, path) in kwargs.iteritems():
            if key in self._paths:
                raise ParameterError("Attempted to overwrite existing path: %s" % key)
        self._paths.update(kwargs)       


    @property
    def call(self):
        """Returns the system-call, based on the call passed to the constructor, and
        every parameter set or pushed using 'set_parameter' and 'push_parameter."""

        command = list(self._call)
        for parameter in self._params:
            if parameter["Value"] is not None:
                if parameter["Sep"] is not None:
                    command.append("%s%s%s" % (parameter["Key"], parameter["Sep"], parameter["Value"]))
                else:
                    command.append(parameter["Key"])
                    command.append(parameter["Value"])
            else:
                command.append(parameter["Key"])
        
        return command


    @property
    def paths(self):
        """Returns a dictionary of paths as set by 'set_paths'."""

        paths = {}
        for (key, value) in self._paths.iteritems():
            if isinstance(value, AtomicParams):
                value = value.finalize()
            paths[key] = value
        return paths


    def finalize(self):
        """Creates an AtomicCmd object based on the AtomicParam object."""
        if not self._cmd_object:
            self._cmd_object = AtomicCmd(self.call, **self.paths)

        return self._cmd_object



class AtomicJavaParams(AtomicParams):
    def __init__(self, config, jar, gc_threads = 1):
        call = ["java", "-server", "-Xmx4g", 
                "-Djava.io.tmpdir=%s" % config.temp_root]

        if gc_threads > 1:
            call.append("-XX:ParallelGCThreads=%i" % gc_threads)
        else:
            call.append("-XX:+UseSerialGC")
            
        call.extend(("-jar", jar))
        AtomicParams.__init__(self, call)



def use_customizable_cli_parameters(init_func):
    """Decorator implementing the customizable Node interface.
    Allows a node to be implemented either using default behavior:
      >>> node = SomeNode(value1 = ..., value2 = ...)
      
    Or using tweaked parameters for calls that support it:
      >>> parameters = SomeNode.customize(value1 = ..., value2 = ...)
      >>> parameters["command"].set_parameter(...)
      >>> node = SomeNode(parameters)

    To be able to use this interface, the class must implement a 
    function 'customize' that takes the parameters that the constructor
    would take, while the constructor must take a single 'parameters'
    key-argument defaulting to None."""

    def do_call(self, parameters = None, **kwargs):
        if not parameters:
            parameters = self.customize(**kwargs)

        return init_func(self, parameters)
    return do_call



def create_customizable_cli_parameters(customize_func):
    def do_call(cls, **kwargs):
        kwargs.update(customize_func(cls, **kwargs))

        # Ensure that parameters with default arguments (e.g. 'dependencies')
        # are passed, if not explicity specified
        spec = inspect.getargspec(customize_func)
        for (key, value) in itertools.izip_longest(reversed(spec.args), reversed(spec.defaults)):
            if key not in kwargs:
                kwargs[key] = value

        name = "%s_parameters" % cls.__name__
        clsobj = collections.namedtuple(name, " ".join(kwargs))
        class _ParametersWrapper(clsobj):
            def build_node(self):
                return cls(self)
        return _ParametersWrapper(**kwargs)
    
    return classmethod(do_call)


