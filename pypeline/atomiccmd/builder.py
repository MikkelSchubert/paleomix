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
        # Use passed parameters to create AtomicCmdBuilder obj
        builder = AtomicCmdBuilder(...)
        builder.set_option(...)

        # Return dictionary of AtomicCmdBuilder objects and any
        # additional parameters required to run the Node.
        return {"command"   : builder,
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

2) Retrieving and tweaking AtomicCmdBuilder before creating the Node:
>> params = ExampleNode.customize(...)
>> params.command.set_option(...)
>> node = params.build_node()

"""
import types
import inspect
import itertools
import collections

from pypeline.atomiccmd.command import AtomicCmd
from pypeline.common.utilities import safe_coerce_to_tuple




class AtomicCmdBuilderError(RuntimeError):
    pass




class AtomicCmdBuilder:
    """AtomicCmdBuilder is a class used to allow step-wise construction of an
    AtomicCmd object. This allows the user of a Node to modify the behavior
    of the called programs using some CLI parameters, without explicit support
    for these in the Node API. Some limitations are in place, to help catch
    cases where overwriting or adding a flag would break the Node.

    The system call is constructed in the following manner:
       $ <call> <options> <values>

    The components are defined as follows:
      <call>    - The minimal call needed invoke the current program. Typically
                  this is just the name of the executable, but may be a more
                  complex set of values for nested calls (e.g. java/scripts).
      <option>  - A flag, typically prefixed with one or two dashes, followed
                  by an optional value. The flag may be joined with the value
                  by an seperator (e.g. '='), otherwise they are added to the
                  final call as seperate values.
      <values>  - One or more values, e.g. paths or similar.

    Options are divided into two classes; singletons and non-singletons:
      Singletons     - May be specified exactly one (using 'set_option'), with
                       subsequent calls to 'set_option' overwriting the previous
                       value of the option (if any).
      Non-singletons - May be specified one or more times (using 'push_option'),
                       with each subsequent call added to list of existing
                       parameters.

    Furthermore, any parameter may be marked as fixed (typically because the
    node dependends on that option being set), which prevents any subsequent
    call from modifying this option. By default, all options are fixed.

    Any number of keywords may be set, which are passed to the AtomicCmd object
    created by the AtomicCmdBuilder object (using 'set_kwargs'). The rules specified
    in the AtomicCmd documentation apply to these. If a AtomicCmdBuilder object
    is passed, this will be finalized as well."""

    def __init__(self, call, **kwargs):
        """See AtomiCmd.__init__ for parameters / keyword arguments."""

        self._call    = safe_coerce_to_tuple(call)
        self._options = []
        self._values  = []
        self._kwargs  = {}
        self._object  = None

        self.set_kwargs(**kwargs)


    def set_option(self, key, value = None, sep = None, fixed = True):
        """Sets or overwrites an option that may be specified at most once. If
        the option has already been set using 'add_option', or with 'fixed' set
        to True, a AtomicCmdBuilderError will be raised."""
        old_option = self._get_option_for_editing(key, singleton = True)
        new_option = {"Key" : key, "Value"  : value, "Sep" : sep, "Fixed"  : fixed, "Singleton" : True}
        if old_option:
            if old_option["Fixed"]:
                raise AtomicCmdBuilderError("Attemping to overwrite fixed option: %r" % key)
            old_option.update(new_option)
        else:
            self._options.append(new_option)


    def add_option(self, key, value = None, sep = None, fixed = True):
        """Adds an option that may be specified one or more times. If the option
        has already been set using 'set_option', a AtomicCmdBuilderError will be raised."""
        # Previous values are not used, but checks are required
        self._get_option_for_editing(key, singleton = False)
        self._options.append({"Key" : key, "Value"  : value, "Sep" : sep, "Fixed"  : fixed, "Singleton" : False})


    def pop_option(self, key):
        old_option = self._get_option_for_editing(key, singleton = None)
        if not old_option:
            raise KeyError("Option with key %r does not exist" % key)
        elif old_option["Fixed"]:
            raise AtomicCmdBuilderError("Attempting to pop fixed key %r" % key)
        self._options.remove(old_option)


    def add_value(self, value):
        """Adds a positional value to the call. Usage should be restricted to
        paths and similar values, and set/add_option used for actual options."""
        self._values.append(value)


    def set_kwargs(self, **kwargs):
        if self._object:
            raise AtomicCmdBuilderError("Parameters have already been finalized")

        for key in kwargs:
            if key in self._kwargs:
                raise AtomicCmdBuilderError("Attempted to overwrite existing path: %r" % key)
        self._kwargs.update(kwargs)


    @property
    def call(self):
        """Returns the system-call, based on the call passed to the constructor, and
        every parameter set or pushed using 'set_option' and 'add_option."""
        command = list(self._call)
        for parameter in self._options:
            if parameter["Value"] is not None:
                if parameter["Sep"] is not None:
                    command.append("%s%s%s" % (parameter["Key"], parameter["Sep"], parameter["Value"]))
                else:
                    command.append(parameter["Key"])
                    command.append(parameter["Value"])
            else:
                command.append(parameter["Key"])

        command.extend(self._values)
        return command


    @property
    def finalized_call(self):
        """Returns the system-call, as 'call', but with all key-values instantiated
        to the values passed to the AtomicCmdBuilder. This is intended for use with
        direct Popen calls."""
        kwargs = self.kwargs
        return [(field % kwargs) for field in self.call]


    @property
    def kwargs(self):
        """Returns a dictionary of keyword arguments as set by 'set_kwargs'.
        If the value of an argument is an AtomicCmdBuilder, then the builder
        is finalized and the resulting value is used."""
        kwargs = {}
        for (key, value) in self._kwargs.iteritems():
            if isinstance(value, AtomicCmdBuilder):
                value = value.finalize()
            kwargs[key] = value
        return kwargs


    def finalize(self):
        """Creates an AtomicCmd object based on the AtomicParam object. Once
        finalized, the AtomicCmdBuilder cannot be modified further."""
        if not self._object:
            self._object = AtomicCmd(self.call, **self.kwargs)

        return self._object


    def _get_option_for_editing(self, key, singleton):
        if self._object:
            raise AtomicCmdBuilderError("AtomicCmdBuilder has already been finalized")
        elif not isinstance(key, types.StringTypes):
            raise TypeError("Key must be a string, not %r" % (key.__class__.__name__,))
        elif not key:
            raise KeyError("Key cannot be an empty string")

        for option in reversed(self._options):
            if (option["Key"] == key):
                if (singleton is not None) and (option["Singleton"] != singleton):
                    raise AtomicCmdBuilderError("Mixing of singleton and non-singleton options: %r" % key)
                return option




class AtomicJavaCmdBuilder(AtomicCmdBuilder):
    def __init__(self, config, jar, gc_threads = 1, **kwargs):
        temp_root = "%(TEMP_DIR)s"
        if config:
            temp_root = config.temp_root

        call = ["java", "-server", "-Xmx4g",
                "-Djava.io.tmpdir=%s" % temp_root,
                "-Djava.awt.headless=true"]

        if not isinstance(gc_threads, (types.IntType, types.LongType)):
            raise TypeError("'gc_threads' must be an integer value, not %r" % gc_threads.__class__.__name__)
        elif gc_threads > 1:
            call.append("-XX:ParallelGCThreads=%i" % gc_threads)
        elif gc_threads == 1:
            call.append("-XX:+UseSerialGC")
        else:
            raise ValueError("'gc_threads' must be a 1 or greater, not %r" % gc_threads)

        call.extend(("-jar", "%(AUX_JAR)s"))
        AtomicCmdBuilder.__init__(self, call, AUX_JAR = jar, **kwargs)




class AtomicMPICmdBuilder(AtomicCmdBuilder):
    def __init__(self, call, threads = 1, **kwargs):
        if not isinstance(threads, (types.IntType, types.LongType)):
            raise TypeError("'threads' must be an integer value, not %r" % threads.__class__.__name__)
        elif threads < 1:
            raise ValueError("'threads' must be 1 or greater, not %i" % threads)
        elif threads == 1:
            AtomicCmdBuilder.__init__(self, call, **kwargs)
        else:
            call = safe_coerce_to_tuple(call)
            mpi_call = ["mpirun", "-n", threads]
            mpi_call.extend(call)

            AtomicCmdBuilder.__init__(self, mpi_call, EXEC_MAIN = call[0], **kwargs)




def use_customizable_cli_parameters(init_func): # pylint: disable=C0103
    """Decorator implementing the customizable Node interface.
    Allows a node to be implemented either using default behavior:
      >>> node = SomeNode(value1 = ..., value2 = ...)

    Or using tweaked parameters for calls that support it:
      >>> parameters = SomeNode.customize(value1 = ..., value2 = ...)
      >>> parameters["command"].set_options(...)
      >>> node = SomeNode(parameters)

    To be able to use this interface, the class must implement a
    function 'customize' that takes the parameters that the constructor
    would take, while the constructor must take a  'parameters' argument."""
    def do_call(self, parameters = None, **kwargs):
        if not parameters:
            parameters = self.customize(**kwargs)

        return init_func(self, parameters)
    return do_call


def create_customizable_cli_parameters(customize_func): # pylint: disable=C0103
    # Ensure that parameters with default arguments (e.g. 'dependencies')
    # are passed, if not explicity specified by the caller
    spec     = inspect.getargspec(customize_func)
    args     = list(reversed(spec.args))
    defaults = list(reversed(spec.defaults or ()))

    def do_call(cls, **kwargs):
        # Allow parameters to be updated in the 'customize' function
        kwargs.update(customize_func(cls, **kwargs))

        for (key, value) in itertools.izip_longest(args, defaults):
            if key not in kwargs:
                kwargs[key] = value

        clsobj = _create_cli_parameters_cls(cls, kwargs)
        return clsobj(**kwargs)

    return classmethod(do_call)



def apply_options(builder, options, pred = lambda s: s.startswith("-")):
    """Applies a dictionary of options to a builder. By default, only
    options where the key start with "-" are used (determined by 'pred').
    The following rules are used when applying options:
      - If a key is assosiated with a single value, 'set_option' is used.
      - If a key is assosiated with a list of values, 'add_option' is used.
      - If the key is assosiated with a boolean value, the option is set
        if true (without a value) or removed from the call if false. This
        allows easy setting/unsetting of '--do-something' type options."""
    for (key, values) in dict(options).iteritems():
        if not isinstance(key, types.StringTypes):
            raise TypeError("Keys must be strings, not %r" % (key.__class__.__name__,))
        elif pred(key):
            if isinstance(values, (types.ListType, types.TupleType)):
                for value in values:
                    if not isinstance(value, _ADDABLE_TYPES) or isinstance(value, _SETABLE_ONLY_TYPES):
                        raise TypeError("Unexpected type when adding options: %r" % (value.__class__.__name__,))
                    builder.add_option(key, value)
            elif not isinstance(values, _SETABLE_TYPES):
                raise TypeError("Unexpected type when setting option: %r" % (values.__class__.__name__,))
            elif isinstance(values, (types.BooleanType, types.NoneType)):
                if values or values is None:
                    builder.set_option(key)
                else:
                    builder.pop_option(key)
            else:
                builder.set_option(key, values)


_create_cli_parameters_cls_cache = {}
def _create_cli_parameters_cls(cls, kwargs):
    key    = (cls, frozenset(kwargs))
    clsobj = _create_cli_parameters_cls_cache.get(key)
    if not clsobj:
        _create_cli_parameters_cls_cache[key] = clsobj = \
          collections.namedtuple("CustomCLIParams", " ".join(kwargs))

    class _ParametersWrapper(clsobj): # pylint: disable=W0232
        def build_node(self):
            return cls(self)

    return _ParametersWrapper



_ADDABLE_TYPES = (types.FloatType, types.IntType, types.LongType) + types.StringTypes
_SETABLE_ONLY_TYPES = (types.BooleanType, types.NoneType)
_SETABLE_TYPES = _ADDABLE_TYPES + _SETABLE_ONLY_TYPES
