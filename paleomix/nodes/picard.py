#!/usr/bin/env python3
"""
Picard - A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.

https://broadinstitute.github.io/picard/
"""
import os
import getpass

import paleomix.common.versions as versions
import paleomix.common.system

from paleomix.atomiccmd.command import AtomicCmd, AuxilleryFile, InputFile, OutputFile
from paleomix.node import CommandNode
from paleomix.common.fileutils import swap_ext, try_rmtree


class PicardNode(CommandNode):
    """Base class for nodes using Picard Tools; adds an additional cleanup
    step, in order to allow the jars to be run using the same temporary folder
    as any other commands associated with the node. This is nessesary as some
    Picard tools create a large number of temporary files, leading to potential
    performance issues if these are located in the same folder.
    """

    def _teardown(self, config, temp):
        # Picard creates a folder named after the user in the temp-root
        try_rmtree(os.path.join(temp, getpass.getuser()))
        # Some JREs may create a folder for temporary performance counters
        try_rmtree(os.path.join(temp, "hsperfdata_" + getpass.getuser()))

        CommandNode._teardown(self, config, temp)


class ValidateBAMNode(PicardNode):
    def __init__(
        self,
        config,
        input_bam,
        input_index=None,
        output_log=None,
        ignored_checks=(),
        big_genome_mode=False,
        dependencies=(),
    ):
        command = _picard_command(
            config,
            [
                "ValidateSamFile",
                "--INPUT",
                InputFile(input_bam),
                "--OUTPUT",
                OutputFile(output_log or swap_ext(input_bam, ".validated")),
            ],
        )

        for check in ignored_checks:
            command.append("--IGNORE", check)

        if big_genome_mode:
            self._configure_for_big_genome(config, command)

        _set_max_open_files(command)

        if input_index is not None:
            command.add_extra_files([InputFile(input_index)])

        PicardNode.__init__(
            self,
            command=command,
            description="validating %s" % (input_bam,),
            dependencies=dependencies,
        )

    @staticmethod
    def _configure_for_big_genome(config, command):
        # CSI uses a different method for assigning BINs to records, which
        # Picard currently does not support.
        command.append("--IGNORE", "INVALID_INDEXING_BIN")

        jar_path = os.path.join(config.jar_root, _PICARD_JAR)
        version_check = _PICARD_VERSION_CACHE[jar_path]

        try:
            if version_check.version >= (2, 19, 0):
                # Useless warning, as we do not build BAI indexes for large genomes
                command.append("--IGNORE", "REF_SEQ_TOO_LONG_FOR_BAI")
        except versions.VersionRequirementError:
            pass  # Ignored here, handled elsewhere


###############################################################################

_PICARD_JAR = "picard.jar"
_PICARD_VERSION_CACHE = {}


def _picard_command(config, args):
    """Returns basic AtomicJavaCmdBuilder for Picard tools commands."""
    jar_path = os.path.join(config.jar_root, _PICARD_JAR)

    if jar_path not in _PICARD_VERSION_CACHE:
        command = _java_cmd(
            jar_path,
            temp_root=config.temp_root,
            jre_options=config.jre_options,
        )

        # Arbitrary command, since just '--version' does not work
        command.append("MarkDuplicates", "--version")

        requirement = versions.Requirement(
            call=command.to_call("%(TEMP_DIR)s"),
            name="Picard tools",
            search=r"\b(\d+)\.(\d+)\.(\d+)",
            checks=versions.GE(2, 10, 8),
        )
        _PICARD_VERSION_CACHE[jar_path] = requirement

    version = _PICARD_VERSION_CACHE[jar_path]
    command = _java_cmd(
        jar_path,
        temp_root=config.temp_root,
        jre_options=config.jre_options,
        requirements=[version],
        set_cwd=True,
    )

    command.append(*args)

    return command


def _java_cmd(jar, jre_options=(), temp_root="%(TEMP_DIR)s", gc_threads=1, **kwargs):
    call = [
        "java",
        "-server",
        "-Djava.io.tmpdir=%s" % temp_root,
        "-Djava.awt.headless=true",
        "-Dpicard.useLegacyParser=false",
    ]

    if not isinstance(gc_threads, int):
        raise TypeError("'gc_threads' must be an integer value, not %r" % (gc_threads,))
    elif gc_threads > 1:
        call.append("-XX:ParallelGCThreads=%i" % gc_threads)
    elif gc_threads == 1:
        call.append("-XX:+UseSerialGC")
    else:
        raise ValueError("'gc_threads' must be a 1 or greater, not %r" % gc_threads)

    call.extend(jre_options)

    # Only set -Xmx if no user-supplied setting is given
    if not any(opt.startswith("-Xmx") for opt in call):
        # Our experience is that the default -Xmx value tends to cause OutOfMemory
        # exceptions with typical datasets, so require at least 4gb.
        call.append("-Xmx4g")

    call.extend(("-jar", AuxilleryFile(jar)))

    return AtomicCmd(call, **kwargs)


# Fraction of per-process max open files to use
_FRAC_MAX_OPEN_FILES = 0.95
# Default maximum number of open temporary files used by Picard
_DEFAULT_MAX_OPEN_FILES = 8000


def _set_max_open_files(command):
    """Sets the maximum number of open files a picard process
    should use, at most. Conservatively lowered than the actual
    ulimit.
    """
    max_open_files = paleomix.common.system.get_max_open_files()
    if max_open_files:
        max_open_files = int(max_open_files * _FRAC_MAX_OPEN_FILES)

        if max_open_files < _DEFAULT_MAX_OPEN_FILES:
            command.append("--MAX_OPEN_TEMP_FILES", max_open_files)
