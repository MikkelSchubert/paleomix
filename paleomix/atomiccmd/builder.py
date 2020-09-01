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
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.common.utilities import safe_coerce_to_tuple


class AtomicCmdBuilderError(RuntimeError):
    """Error raised by AtomicCmdBuilder."""


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
                       subsequent calls to 'set_option' overwriting the
                       previous value of the option (if any).
      Non-singletons - May be specified one or more times (with 'push_option'),
                       with each subsequent call added to list of existing
                       parameters.

    Furthermore, any parameter may be marked as fixed (typically because the
    node dependends on that option being set), which prevents any subsequent
    call from modifying this option. By default, all options are fixed.

    Any number of keywords may be set, which are passed to the AtomicCmd object
    created by the AtomicCmdBuilder object (using 'set_kwargs'). The rules
    specified in the AtomicCmd documentation apply to these. If a
    AtomicCmdBuilder object is passed, this will be finalized as well.
    """

    def __init__(self, call, **kwargs):
        """See AtomiCmd.__init__ for parameters / keyword arguments.
        """
        self._call = safe_coerce_to_tuple(call)
        self._options = []
        self._values = []
        self._kwargs = {}
        self._object = None

        self.set_kwargs(**kwargs)

    def set_option(self, key, value=None, sep=None, fixed=True):
        """Sets or overwrites an option that may be specified at most once. If
        the option has already been set using 'add_option', or with 'fixed' set
        to True, a AtomicCmdBuilderError will be raised."""
        old_option = self._get_option_for_editing(key, singleton=True)
        new_option = {
            "Key": key,
            "Value": value,
            "Sep": sep,
            "Fixed": fixed,
            "Singleton": True,
        }

        if old_option:
            if old_option["Fixed"]:
                message = "Attemping to overwrite fixed option: %r" % key
                raise AtomicCmdBuilderError(message)
            old_option.update(new_option)
        else:
            self._options.append(new_option)

    def add_option(self, key, value=None, sep=None, fixed=True):
        """Adds an option that may be specified one or more times. If the
        option has already been set using 'set_option', a AtomicCmdBuilderError
        will be raised.
        """
        # Previous values are not used, but checks are required
        self._get_option_for_editing(key, singleton=False)
        self._options.append(
            {"Key": key, "Value": value, "Sep": sep, "Fixed": fixed, "Singleton": False}
        )

    def pop_option(self, key):
        """Remove option with key; raises error if option does not exist."""
        old_option = self._get_option_for_editing(key, singleton=None)
        if not old_option:
            raise KeyError("Option with key %r does not exist" % key)
        elif old_option["Fixed"]:
            raise AtomicCmdBuilderError("Attempting to pop fixed key %r" % key)
        self._options.remove(old_option)

    def add_value(self, value):
        """Adds a positional value to the call. Usage should be restricted to
        paths and similar values, and set/add_option used for actual options.
        """
        self._values.append(value)

    def set_kwargs(self, **kwargs):
        """Sets any number of keyword arguments; raises an exception if any of
        the arguments have already been set."""
        if self._object:
            message = "Parameters have already been finalized"
            raise AtomicCmdBuilderError(message)

        for key in kwargs:
            if key in self._kwargs:
                message = "Attempted to overwrite existing path: %r"
                raise AtomicCmdBuilderError(message % key)
        self._kwargs.update(kwargs)

    def add_multiple_options(self, key, values, sep=None, template="IN_FILE_%02i"):
        """Add multiple options as once, with corresponding kwargs.

        The template determines the key-names used for the arguments,
        using numbers starting from 1 to differentiate between multiple
        values.
        """
        kwargs = {}
        for file_key, value in self._get_new_kwarg_keys(values, template):
            self.add_option(key, "%%(%s)s" % (file_key,), sep=sep, fixed=True)
            kwargs[file_key] = value
        self.set_kwargs(**kwargs)
        return kwargs

    def add_multiple_values(self, values, template="IN_FILE_%02i"):
        """Add multiple values as once, with corresponding kwargs.

        The template determines the key-names used for the arguments,
        using numbers starting from 1 to differentiate between multiple
        values.
        """
        kwargs = {}
        for file_key, value in self._get_new_kwarg_keys(values, template):
            self.add_value("%%(%s)s" % (file_key,))
            kwargs[file_key] = value
        self.set_kwargs(**kwargs)
        return kwargs

    def add_multiple_kwargs(self, values, template="IN_FILE_%02i"):
        """Add multiple keyword arguments using the given values.

        The template determines the key-names used for the arguments,
        using numbers starting from 1 to differentiate between multiple
        values.
        """
        kwargs = dict(self._get_new_kwarg_keys(values, template))
        self.set_kwargs(**kwargs)
        return kwargs

    @property
    def call(self):
        """Returns the system-call based on the call passed to the constructor,
        and every parameter set or pushed using 'set_option' and 'add_option'.
        """
        command = list(self._call)
        for parameter in self._options:
            if parameter["Value"] is not None:
                if parameter["Sep"] is not None:
                    command.append(
                        "%s%s%s"
                        % (parameter["Key"], parameter["Sep"], parameter["Value"])
                    )
                else:
                    command.append(parameter["Key"])
                    command.append(parameter["Value"])
            else:
                command.append(parameter["Key"])

        command.extend(self._values)
        return command

    @property
    def finalized_call(self):
        """Returns the system-call, as 'call', but with all key-values
        instantiated to the values passed to the AtomicCmdBuilder. This is
        intended for use with direct Popen calls.
        """
        kwargs = self.kwargs
        kwargs["TEMP_DIR"] = "%(TEMP_DIR)"
        return [(str(field) % kwargs) for field in self.call]

    @property
    def kwargs(self):
        """Returns a dictionary of keyword arguments as set by 'set_kwargs'.
        If the value of an argument is an AtomicCmdBuilder, then the builder
        is finalized and the resulting value is used."""
        kwargs = {}
        for (key, value) in self._kwargs.items():
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
            message = "AtomicCmdBuilder has already been finalized"
            raise AtomicCmdBuilderError(message)
        elif not isinstance(key, str):
            message = "Key must be a string, not %r" % (key.__class__.__name__,)
            raise TypeError(message)
        elif not key:
            raise KeyError("Key cannot be an empty string")

        for option in reversed(self._options):
            if option["Key"] == key:
                if (singleton is not None) and (option["Singleton"] != singleton):
                    message = "Mixing singleton and non-singleton options: %r"
                    raise AtomicCmdBuilderError(message % key)
                return option

    def _get_new_kwarg_keys(self, values, template):
        start = 0
        for value in values:
            start += 1
            key = template % (start,)
            while key in self._kwargs:
                start += 1
                key = template % (start,)
            yield key, value


class AtomicJavaCmdBuilder(AtomicCmdBuilder):
    """AtomicCmdBuilder for running java JARs.

    The resulting command will run the JAR in head-less mode, in order to
    ensure that the JARs can be run on head-less servers (and to avoid popups
    on OSX), using the process-specific temp-folder, and using at most a single
    thread for garbage collection (to ensure that thread-limits are obeyed).
    """

    def __init__(
        self, jar, jre_options=(), temp_root="%(TEMP_DIR)s", gc_threads=1, **kwargs
    ):
        """Parameters:
            jar         -- Path to a JAR file to be executed; is included as an
                           auxiliary file dependency in the final command.
            jre_options -- List of CLI options to be passed to 'java' command.
            temp_root   -- Temp folder to use for java process; if not set, the
                           process specific temp folder is used.
            gc_threads  -- Number of threads to use during garbage collections.
            ...         -- Key-word args are passed to AtomicCmdBuilder.
        """
        call = [
            "java",
            "-server",
            "-Djava.io.tmpdir=%s" % temp_root,
            "-Djava.awt.headless=true",
        ]

        if not isinstance(gc_threads, int):
            raise TypeError(
                "'gc_threads' must be an integer value, not %r"
                % gc_threads.__class__.__name__
            )
        elif gc_threads > 1:
            call.append("-XX:ParallelGCThreads=%i" % gc_threads)
        elif gc_threads == 1:
            call.append("-XX:+UseSerialGC")
        else:
            raise ValueError("'gc_threads' must be a 1 or greater, not %r" % gc_threads)

        jre_options = tuple(jre_options)
        call.extend(jre_options)

        # Only set -Xmx if no user-supplied setting is given
        if not any(opt.startswith("-Xmx") for opt in jre_options):
            # Our experience is that the default -Xmx value tends to cause OutOfMemory
            # exceptions with typical datasets, so require at least 4gb.
            call.append("-Xmx4g")

        call.extend(("-jar", "%(AUX_JAR)s"))
        AtomicCmdBuilder.__init__(self, call, AUX_JAR=jar, **kwargs)


class AtomicMPICmdBuilder(AtomicCmdBuilder):
    """AtomicCmdBuilder for MPI enabled programs;

    Simplifies specification of number of threads to use, only invoking the
    'mpi' command if more than one thread is used; furthermore, the 'mpi'
    binary is used as a dependency, since MPI enabled programs tend to fail
    catastrophically if the 'mpi' binary and associated libraries are missing.

    """

    def __init__(self, call, threads=1, **kwargs):
        if not isinstance(threads, int):
            raise TypeError(
                "'threads' must be an integer value, not %r"
                % threads.__class__.__name__
            )
        elif threads < 1:
            raise ValueError("'threads' must be 1 or greater, not %i" % threads)
        elif threads == 1:
            AtomicCmdBuilder.__init__(self, call, EXEC_MPI="mpirun", **kwargs)
        else:
            call = safe_coerce_to_tuple(call)
            mpi_call = ["mpirun", "-n", threads]
            mpi_call.extend(call)

            AtomicCmdBuilder.__init__(self, mpi_call, EXEC_MAIN=call[0], **kwargs)


def apply_options(builder, options, pred=lambda s: s.startswith("-")):
    """Applies a dictionary of options to a builder. By default, only
    options where the key start with "-" are used (determined by 'pred').
    The following rules are used when applying options:
      - If a key is associated with a single value, 'set_option' is used.
      - If a key is associated with a list of values, 'add_option' is used.
      - If the key is associated with a boolean value, the option is set
        if true (without a value) or removed from the call if false. This
        allows easy setting/unsetting of '--do-something' type options.

    """
    for (key, values) in dict(options).items():
        if not isinstance(key, str):
            raise TypeError("Keys must be strings, not %r" % (key.__class__.__name__,))
        elif pred(key):
            if isinstance(values, (list, tuple)):
                for value in values:
                    if not isinstance(value, _ADDABLE_TYPES) or isinstance(
                        value, _SETABLE_ONLY_TYPES
                    ):
                        raise TypeError(
                            "Unexpected type when in options: %r"
                            % (value.__class__.__name__,)
                        )
                    builder.add_option(key, value)
            elif not isinstance(values, _SETABLE_TYPES):
                raise TypeError(
                    "Unexpected type when setting option: %r"
                    % (values.__class__.__name__,)
                )
            elif isinstance(values, (bool, type(None))):
                if values or values is None:
                    builder.set_option(key)
                else:
                    builder.pop_option(key)
            else:
                builder.set_option(key, values)


_ADDABLE_TYPES = (float, int, str)
_SETABLE_ONLY_TYPES = (bool, type(None))
_SETABLE_TYPES = _ADDABLE_TYPES + _SETABLE_ONLY_TYPES
