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
import copy

from paleomix.node import Node
from paleomix.common.fileutils import move_file, reroot_path, describe_files
from paleomix.common.formats.msa import MSA
from paleomix.common.formats.phylip import interleaved_phy

from paleomix.common.utilities import safe_coerce_to_frozenset, safe_coerce_to_tuple


_VALID_KEYS = frozenset(["partitions", "filenames"])


class FastaToPartitionedInterleavedPhyNode(Node):
    def __init__(
        self,
        infiles,
        out_prefix,
        exclude_groups=(),
        reduce=False,
        dependencies=(),
        file_dependencies=(),
    ):
        """
        infiles = {names : {"partitions" : ..., "filenames" : [...]}}
        """
        if not (
            isinstance(infiles, dict)
            and all(isinstance(dd, dict) for dd in infiles.values())
        ):
            raise TypeError("'infiles' must be a dictionary of dictionaries")

        input_filenames = []
        for (name, subdd) in infiles.items():
            if set(subdd) - _VALID_KEYS:
                raise ValueError(
                    "Invalid keys found for %r: %s"
                    % (name, ", ".join(set(subdd) - _VALID_KEYS))
                )
            elif not isinstance(subdd["filenames"], list):
                raise ValueError("filenames must be a list of strings")
            input_filenames.extend(subdd["filenames"])
        # Optional file dependencies; used to depend on the list of sequcences
        input_filenames.extend(safe_coerce_to_tuple(file_dependencies))

        self._reduce = bool(reduce)
        self._infiles = copy.deepcopy(infiles)
        self._out_prefix = out_prefix
        self._excluded = safe_coerce_to_frozenset(exclude_groups)

        description = "creating%spartitioned phy from %s" % (
            "  reduced, " if reduce else " ",
            describe_files(input_filenames),
        )

        Node.__init__(
            self,
            description=description,
            input_files=input_filenames,
            output_files=[out_prefix + ".phy", out_prefix + ".partitions"],
            dependencies=dependencies,
        )

    def _run(self, _config, temp):
        merged_msas = []
        for (name, files_dd) in sorted(self._infiles.items()):
            partitions = files_dd["partitions"]
            msas = dict((key, []) for key in partitions)
            for filename in files_dd["filenames"]:
                msa = MSA.from_file(filename)
                if self._excluded:
                    msa = msa.exclude(self._excluded)

                for (key, msa_part) in msa.split(partitions).items():
                    msas[key].append(msa_part)

            msas.pop("X", None)
            for (key, msa_parts) in sorted(msas.items()):
                merged_msa = MSA.join(*msa_parts)
                if self._reduce:
                    merged_msa = merged_msa.reduce()

                if merged_msa is not None:
                    merged_msas.append(("%s_%s" % (name, key), merged_msa))

        out_fname_phy = reroot_path(temp, self._out_prefix + ".phy")
        with open(out_fname_phy, "w") as output_phy:
            final_msa = MSA.join(*(msa for (_, msa) in merged_msas))
            output_phy.write(interleaved_phy(final_msa))

        partition_end = 0
        out_fname_parts = reroot_path(temp, self._out_prefix + ".partitions")
        with open(out_fname_parts, "w") as output_part:
            for (name, msa) in merged_msas:
                length = msa.seqlen()
                output_part.write(
                    "DNA, %s = %i-%i\n"
                    % (name, partition_end + 1, partition_end + length)
                )
                partition_end += length

    def _teardown(self, _config, temp):
        move_file(
            reroot_path(temp, self._out_prefix + ".phy"), self._out_prefix + ".phy"
        )
        move_file(
            reroot_path(temp, self._out_prefix + ".partitions"),
            self._out_prefix + ".partitions",
        )
