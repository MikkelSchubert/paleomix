#
# Copyright (c) 2023 Mikkel Schubert <MikkelSch@gmail.com>
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
from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from paleomix.common import fileutils, versions
from paleomix.common.command import (
    AtomicCmd,
    InputFile,
    OptionsType,
    OutputFile,
    TempOutputFile,
)
from paleomix.node import CommandNode, Node

_VERSION_2_CHECK = versions.Requirement(
    call=("AdapterRemoval", "--version"),
    regexp=r"ver. (\d+\.\d+\.\d+)",
    specifiers=">=2.2.0",
)


_VERSION_3_CHECK = versions.Requirement(
    call=("adapterremoval3", "--version"),
    # Currently fixed to specify development version
    regexp=r"AdapterRemoval v(\d+\.\d+\.\d+-alpha\d+)",
    specifiers=">=3.0.0-alpha2",
)


class AdapterRemoval2Node(CommandNode):
    out_fastq: dict[str, tuple[str, str | None]]

    def __init__(
        self,
        input_file_1: str,
        input_file_2: str | None,
        output_prefix: str,
        output_settings: str,
        threads: int = 1,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        # Options that cannot be overwritten
        fixed_options: OptionsType = {
            "--file1": InputFile(input_file_1),
            # Gzip compress FASTQ files
            "--gzip": None,
            # Fix number of threads to ensure consistency when scheduling node
            "--threads": threads,
            # Prefix for output files, ensure that all end up in temp folder
            "--basename": TempOutputFile(output_prefix),
            # Possibly non-standard locations for settings
            "--settings": OutputFile(output_settings),
        }

        if input_file_2 is None:
            self.out_fastq = {
                "Single": (f"{output_prefix}.truncated.gz", None),
                "Discarded": (f"{output_prefix}.discarded.gz", None),
            }

            # The ability to "collapse" SE reads is not used
            options = {
                key: value
                for key, value in options.items()
                if not key.startswith("--collapse")
            }
        else:
            self.out_fastq = {
                "Paired": (
                    f"{output_prefix}.pair1.truncated.gz",
                    f"{output_prefix}.pair2.truncated.gz",
                ),
                "Singleton": (f"{output_prefix}.singleton.truncated.gz", None),
                "Discarded": (f"{output_prefix}.discarded.gz", None),
            }

            fixed_options["--file2"] = InputFile(input_file_2)

            # merging is enabled if any of the --collapse arguments are used
            if any(key.startswith("--collapse") for key in options):
                self.out_fastq["Collapsed"] = (f"{output_prefix}.collapsed.gz", None)
                self.out_fastq["CollapsedTruncated"] = (
                    f"{output_prefix}.collapsed.truncated.gz",
                    None,
                )

        command = _finalize_options(
            command=AtomicCmd(
                "AdapterRemoval",
                requirements=[_VERSION_2_CHECK],
            ),
            out_fastq=self.out_fastq,
            user_options=options,
            fixed_options=fixed_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            threads=threads,
            description=_describe_trimming_task(input_file_1, input_file_2),
            dependencies=dependencies,
        )


class AdapterRemoval3Node(CommandNode):
    out_fastq: dict[str, tuple[str, str | None]]

    def __init__(
        self,
        input_file_1: str,
        input_file_2: str | None,
        output_prefix: str,
        output_json: str,
        output_html: str,
        threads: int = 1,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        # Options that cannot be overwritten
        fixed_options: OptionsType = {
            "--file1": InputFile(input_file_1),
            # Gzip compress FASTQ files
            "--gzip": None,
            # Fix number of threads to ensure consistency when scheduling node
            "--threads": threads,
            # Prefix for output files, ensure that all end up in temp folder
            "--basename": TempOutputFile(output_prefix),
            # Possibly non-standard locations for settings
            "--out-json": OutputFile(output_json),
            # Possibly non-standard locations for settings
            "--out-html": OutputFile(output_html),
        }

        options = dict(options)
        # Option is no longer supported (output is always Phred+33)
        options.pop("--qualitybase-output", None)

        if input_file_2 is None:
            self.out_fastq = {
                "Single": (f"{output_prefix}.r1.fastq.gz", None),
                "Discarded": (f"{output_prefix}.discarded.fastq.gz", None),
            }
        else:
            self.out_fastq = {
                "Paired": (
                    f"{output_prefix}.r1.fastq.gz",
                    f"{output_prefix}.r2.fastq.gz",
                ),
                "Singleton": (f"{output_prefix}.singleton.fastq.gz", None),
                "Discarded": (f"{output_prefix}.discarded.fastq.gz", None),
            }

            fixed_options["--file2"] = InputFile(input_file_2)

            # merging is enabled if any of the --merge/--collapse arguments are used
            if any(key.startswith(("--collapse", "--merge")) for key in options):
                self.out_fastq["Collapsed"] = (f"{output_prefix}.merged.fastq.gz", None)

        command = _finalize_options(
            command=AtomicCmd(
                "adapterremoval3",
                requirements=[_VERSION_3_CHECK],
            ),
            out_fastq=self.out_fastq,
            user_options=options,
            fixed_options=fixed_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            threads=threads,
            description=_describe_trimming_task(input_file_1, input_file_2),
            dependencies=dependencies,
        )


class IdentifyAdaptersNode(CommandNode):
    def __init__(
        self,
        input_file_1: str,
        input_file_2: str,
        output_file: str,
        threads: int = 1,
        options: OptionsType | None = None,
        dependencies: Iterable[Node] = (),
    ) -> None:
        if options is None:
            options = {}

        command = AtomicCmd(
            "AdapterRemoval",
            stdout=output_file,
            requirements=[_VERSION_2_CHECK],
        )

        fixed_options: OptionsType = {
            "--identify-adapters": None,
            "--file1": InputFile(input_file_1),
            "--file2": InputFile(input_file_2),
            # Fix number of threads to ensure consistency when scheduling node
            "--threads": threads,
        }

        command.merge_options(
            user_options=options,
            fixed_options=fixed_options,
        )

        CommandNode.__init__(
            self,
            command=command,
            threads=threads,
            description="identifying PE adapters in {}".format(
                fileutils.describe_paired_files(input_file_1, input_file_2)
            ),
            dependencies=dependencies,
        )


def _describe_trimming_task(input_file_1: str, input_file_2: str | None) -> str:
    if input_file_2 is None:
        return f"trimming SE adapters from {fileutils.describe_files(input_file_1)}"

    return "trimming PE adapters from {}".format(
        fileutils.describe_paired_files(input_file_1, input_file_2)
    )


def _finalize_options(
    command: AtomicCmd,
    out_fastq: dict[str, tuple[str, str | None]],
    user_options: OptionsType,
    fixed_options: OptionsType,
) -> AtomicCmd:
    output_fastq: list[OutputFile] = []
    for file_1, file_2 in out_fastq.values():
        output_fastq.append(OutputFile(file_1))
        if file_2 is not None:
            output_fastq.append(OutputFile(file_2))

    command.add_extra_files(output_fastq)

    # Ensure that any user-specified list of adapters is tracked
    adapter_list = user_options.get("--adapter-list")
    if isinstance(adapter_list, str):
        user_options["--adapter-list"] = InputFile(adapter_list)

    # Terminal trimming options are specified as --trimXp A B
    extra_options: list[Any] = []
    for key in ("--trim5p", "--trim3p"):
        value = user_options.get(key)
        if isinstance(value, list):
            extra_options.append(key)
            extra_options.extend(value)
            user_options.pop(key)

    command.merge_options(
        user_options=user_options,
        fixed_options=fixed_options,
    )

    command.append(*extra_options)

    return command
