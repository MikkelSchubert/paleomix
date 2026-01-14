# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
"""
BCFtools - BCFtools is a program for variant calling and manipulating files in the
Variant Call Format (VCF) and its binary counterpart BCF.

https://github.com/samtools/bcftools
"""

from __future__ import annotations

from collections.abc import Iterable

from paleomix.common import versions
from paleomix.common.command import (
    AtomicCmd,
    InputFile,
    OptionsType,
    OutputFile,
)
from paleomix.node import CommandNode, Node

_VERSION_REGEX = r"Version: (\d+\.\d+)(?:\.(\d+))?"


# Version >= 1.20 required for --write-index=FMT
BCFTOOLS_VERSION = versions.Requirement(
    call=("bcftools",),
    regexp=_VERSION_REGEX,
    specifiers=">=1.20.0",
)


class BCFMergeNode(CommandNode):
    """Merge VCF or BCF files"""

    def __init__(
        self,
        in_variants: list[str],
        out_variants: str,
        options: OptionsType,
        dependencies: Iterable[Node] = (),
    ) -> None:
        cmd = AtomicCmd(
            ["bcftools", "merge"],
            requirements=[BCFTOOLS_VERSION],
        )

        options = dict(options)
        options["--output"] = OutputFile(out_variants)

        for key in ("-g", "--gvcf"):
            if key in options:
                filename = options[key]
                if isinstance(filename, InputFile):
                    filename = filename.path
                elif not isinstance(filename, str):
                    raise ValueError((key, filename))

                cmd.add_extra_files([InputFile(f"{filename}.fai")])

        index_file: None | OutputFile = None
        for key in options:
            for write_index_key in ("-W", "--write-index"):
                if key.startswith(write_index_key):
                    fmt = key[len(write_index_key) :].lstrip("=")
                    if fmt:
                        index_file = OutputFile(f"{out_variants}.{fmt}")
                    else:
                        raise RuntimeError(key)

        if index_file is not None:
            cmd.add_extra_files([index_file])

        threads = options.setdefault("--threads", 1)
        if not isinstance(threads, int):
            raise ValueError(("--threads", threads))

        cmd.append_options(options)

        for filename in in_variants:
            cmd.append(InputFile(filename))

        CommandNode.__init__(
            self,
            threads=threads,
            description=f"merging files into {out_variants}",
            command=cmd,
            dependencies=dependencies,
        )
