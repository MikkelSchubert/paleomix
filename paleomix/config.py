#!/usr/bin/python
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
import os
import sys
import types
import socket
import getpass
import optparse
import ConfigParser
import multiprocessing

from paleomix.common.fileutils import \
     make_dirs
from paleomix.common.console import \
     print_info


class ConfigError(RuntimeError):
    pass


class PerHostValue:
    def __init__(self, value, is_path=False):
        """Represents a config value that should be settable on a
        per-host basis. If 'is_path' is set, the value is assumed
        to represent a path, and ~ is expanded to the user's home
        folder."""

        self.value = value
        self.is_path = is_path


class PerHostConfig:
    """Helper class for optparse.OptionParser use by pipelines; standardizes the
    process of reading / writing overridable CLI options, while allowing per-host
    options to be set in the .ini files. Values for options with a default value
    that is an instance of PerHostValue will automatically be read-from /
    written-to the per-user config file.

    Given a pipeline with name "NAME", the class will read options from the
    following configuration files:
      - /etc/pypeline/NAME.ini
      - ~/.pypeline/NAME.ini

    These files are expected to contain a "Defaults" section (applies to all hosts),
    and an optional section using the hostname of a server. A file containing the
    current settings (as passed on the CLI) may be written using --write-config-file.

    Example usage:
      per_host_cfg = PerHostConfig("my_pypeline")
      parser       = OptionParser(...)
      parser.add_option(..., default = PerHostValue(...))
      config, args = per_host_cfg.parse_args(parser, sys.argv[1:])
    """

    def __init__(self, pipeline_name):
        """Creates a PerHostConfig for a pipeline with the specified name."""
        # Various common options
        temp_root = os.path.join("/tmp", getpass.getuser(), pipeline_name)
        self.temp_root = PerHostValue(temp_root, True)
        # At least 2 threads are required for e.g. PE BWA nodes, and generally recommended anyway
        self.max_threads = PerHostValue(max(2, multiprocessing.cpu_count()))

        self._filenames = self._get_filenames(pipeline_name)
        self._handle = ConfigParser.SafeConfigParser()
        self._handle.read(self._filenames)
        self._sections = []

        hostname = socket.gethostname()
        if self._handle.has_section(hostname):
            self._sections.append(hostname)
        self._sections.append("Defaults")

    def parse_args(self, parser, argv):
        """Calls 'parse_args' on the parser object after updating default values
        using the settings-files. If --write-config-file is set, a config file
        containing the resulting settings is written."""
        self._add_per_host_options(parser)
        defaults = self._update_defaults(parser)
        config, args = parser.parse_args(argv)

        if config.write_config_file:
            self._write_config_file(config, defaults)

        return config, args

    def _write_config_file(self, config, defaults):
        """Writes a basic config files, using the values previously found in the
        config files, and specified on the command-line."""
        defaults_cfg = ConfigParser.SafeConfigParser()
        defaults_cfg.add_section("Defaults")
        for key in defaults:
            value = getattr(config, key)
            if isinstance(value, (types.ListType, types.TupleType)):
                value = ";".join(value)

            defaults_cfg.set("Defaults", key, str(value))

        filename = self._filenames[-1]
        make_dirs(os.path.dirname(filename))
        with open(filename, "w") as handle:
            defaults_cfg.write(handle)

        print_info("Wrote config file %r" % (filename,))
        sys.exit(0)

    def _add_per_host_options(self, parser):
        """Adds options to a parser relating to the PerHostConfig."""
        group = optparse.OptionGroup(parser, "Config files")
        group.add_option("--write-config-file",
                         default=False, action="store_true",
                         help="Write config using current settings to %s"
                         % (self._filenames[-1],))
        parser.add_option_group(group)

    def _update_defaults(self, parser):
        """Updates default values in a OptionParser, and returns a new
        ConfigParser object containing a new default-values object derived
        from current config-files / CLI options."""
        defaults = {}
        for opt in parser._get_all_options():
            if isinstance(opt.default, PerHostValue):
                defaults[opt.dest] = self._get_default(opt)
        parser.set_defaults(**defaults)
        return defaults

    def _get_default(self, option):
        value = option.default.value
        getter = self._handle.get
        if isinstance(value, types.BooleanType):
            getter = self._handle.getboolean
        elif isinstance(value, (types.IntType, types.LongType)):
            getter = self._handle.getint
        elif isinstance(value, (types.FloatType)):
            getter = self._handle.getfloat
        elif isinstance(value, (types.ListType, types.TupleType)):
            def getter(section, key):
                return filter(None, self._handle.get(section, key).split(";"))

        for section in self._sections:
            if self._handle.has_option(section, option.dest):
                value = getter(section, option.dest)
                break

        if option.default.is_path:
            value = os.path.expanduser(value)

        return value

    @classmethod
    def _get_filenames(cls, name):
        """Return standard list of config files for PALEOMIX pipelines:
           - /etc/pypeline/{name}.ini
           - /etc/paleomix/{name}.ini
           - ~/.paleomix/{name}.ini
        """
        filename = "%s.ini" % (name,)
        homefolder = os.path.expanduser('~')
        return ["/etc/pypeline/%s" % (filename,),
                "/etc/paleomix/%s" % (filename,),
                os.path.join(homefolder, ".paleomix", filename)]


def migrate_config():
    """Checks for the existence of PALEOMIX config files in the old, deprecated
    location (~/.pypeline), and moves these to the new location (~/.paleomix),
    if no config files exist. The old location is replaced with a symbolic
    link, to ensure that older versions of the pipeline do not break.
    """
    homefolder = os.path.expanduser('~')
    old_root = os.path.join(homefolder, ".pypeline")
    new_root = os.path.join(homefolder, ".paleomix")

    if not os.path.exists(new_root):
        if os.path.exists(old_root):
            sys.stderr.write("INFO: Migrating ~/.pypeline to ~/.paleomix\n")
            os.rename(old_root, new_root)
            os.symlink(new_root, old_root)
