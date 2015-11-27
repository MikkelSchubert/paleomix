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
import copy
import collections

from paleomix.node import Node
from paleomix.common.fileutils import move_file, reroot_path
from paleomix.common.formats.msa import MSA
from paleomix.common.formats.phylip import interleaved_phy, sequential_phy

from paleomix.common.utilities import \
     safe_coerce_to_frozenset, \
     safe_coerce_to_tuple



_VALID_KEYS = frozenset(["partitions", "filenames"])


class FastaToPartitionedInterleavedPhyNode(Node):
    def __init__(self, infiles, out_prefix, exclude_groups=(), reduce=False,
                 dependencies=(), file_dependencies=()):
        """
        infiles = {names : {"partitions" : ..., "filenames" : [...]}}
        """
        if not (isinstance(infiles, dict)
                and all(isinstance(dd, dict) for dd in infiles.values())):
            raise TypeError("'infiles' must be a dictionary of dictionaries")

        input_filenames = []
        for (name, subdd) in infiles.iteritems():
            if set(subdd) - _VALID_KEYS:
                raise ValueError("Invalid keys found for %r: %s"
                                 % (name, ", ".join(set(subdd) - _VALID_KEYS)))
            elif not isinstance(subdd["filenames"], list):
                raise ValueError("filenames must be a list of strings")
            input_filenames.extend(subdd["filenames"])
        # Optional file dependencies; used to depend on the list of sequcences
        input_filenames.extend(safe_coerce_to_tuple(file_dependencies))

        self._reduce = bool(reduce)
        self._infiles = copy.deepcopy(infiles)
        self._out_prefix = out_prefix
        self._excluded = safe_coerce_to_frozenset(exclude_groups)

        description = "<FastaToPartitionedPhy%s: %i file(s) -> '%s.*'>" % \
            (" (reducing)" if reduce else "", len(infiles), out_prefix)

        Node.__init__(self,
                      description=description,
                      input_files=input_filenames,
                      output_files=[out_prefix + ".phy",
                                    out_prefix + ".partitions"],
                      dependencies=dependencies)

    def _run(self, _config, temp):
        merged_msas = []
        for (name, files_dd) in sorted(self._infiles.iteritems()):
            partitions = files_dd["partitions"]
            msas = dict((key, []) for key in partitions)
            for filename in files_dd["filenames"]:
                msa = MSA.from_file(filename)
                if self._excluded:
                    msa = msa.exclude(self._excluded)

                for (key, msa_part) in msa.split(partitions).iteritems():
                    msas[key].append(msa_part)

            msas.pop("X", None)
            for (key, msa_parts) in sorted(msas.iteritems()):
                merged_msa = MSA.join(*msa_parts)
                if self._reduce:
                    merged_msa = merged_msa.reduce()

                if merged_msa is not None:
                    merged_msas.append(("%s_%s" % (name, key),
                                        merged_msa))

        out_fname_phy = reroot_path(temp, self._out_prefix + ".phy")
        with open(out_fname_phy, "w") as output_phy:
            final_msa = MSA.join(*(msa for (_, msa) in merged_msas))
            output_phy.write(interleaved_phy(final_msa))

        partition_end = 0
        out_fname_parts = reroot_path(temp, self._out_prefix + ".partitions")
        with open(out_fname_parts, "w") as output_part:
            for (name, msa) in merged_msas:
                length = msa.seqlen()
                output_part.write("DNA, %s = %i-%i\n"
                                  % (name,
                                     partition_end + 1,
                                     partition_end + length))
                partition_end += length

    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_prefix + ".phy"),
                  self._out_prefix + ".phy")
        move_file(reroot_path(temp, self._out_prefix + ".partitions"),
                  self._out_prefix + ".partitions")




class FastaToPartitionsNode(Node):
    def __init__(self, infiles, out_partitions, partition_by = "123", dependencies = ()):
        if (len(partition_by) != 3):
            raise ValueError("Default 'partition_by' must be 3 entires long!")
        elif not isinstance(infiles, dict):
            raise TypeError("'infiles' must be a dictionary")
        elif any(len(dd.get("partition_by", "123")) != 3 for dd in infiles.itervalues()):
            raise ValueError("'partition_by' must be 3 entires long!")
        elif not all(isinstance(dd, dict) for dd in infiles.values()):
            raise TypeError("'infiles' must be a dictionary of dictionaries")
        elif not any(("name" in dd) for dd in infiles.values()):
            raise ValueError("'name' must be specified for all input files")
        elif any((set(dd) - _VALID_KEYS) for dd in infiles.values()):
            raise ValueError("Invalid keys found: %s" % ", ".join(set(dd) - _VALID_KEYS))

        self._infiles   = infiles
        self._out_part  = out_partitions
        self._part_by   = partition_by

        description  = "<FastaToPartitions (default: %s): %i file(s) -> '%s'>" % \
            (partition_by, len(infiles), out_partitions)

        Node.__init__(self,
                      description  = description,
                      input_files  = infiles.keys(),
                      output_files = out_partitions,
                      dependencies = dependencies)


    def _run(self, _config, temp):
        end = 0
        partitions = collections.defaultdict(list)
        for (filename, msa) in _read_sequences(self._infiles):
            length = msa.seqlen()
            start, end = end + 1, end + length

            for (group, offsets) in self._get_partition_by(filename):
                if len(offsets) != 3:
                    parts = [("%i-%i\\3" % (start + offset, end)) for offset in offsets]
                else:
                    parts = ["%i-%i" % (start, end)]

                name = "%s_%s" % (self._infiles[filename]["name"], group)
                partitions[name].extend(parts)

        with open(reroot_path(temp, self._out_part), "w") as part_file:
            for (name, parts) in sorted(partitions.items()):
                part_file.writelines("DNA, %s = %s\n" % (name, ", ".join(parts)))


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_part), self._out_part)


    def _get_partition_by(self, filename):
        groups = self._infiles[filename].get("partition_by", self._part_by)

        partition_by = {}
        for (group, offset) in zip(groups, range(3)):
            partition_by.setdefault(group, []).append(offset)

        return list(sorted(partition_by.items()))



class FastaToInterleavedPhyNode(Node):
    def __init__(self, infiles, out_phy, add_flag = False, dependencies = ()):
        self._add_flag  = add_flag
        self._out_phy   = out_phy

        description  = "<FastaToInterleavedPhy: %i file(s) -> '%s'%s>" % \
            (len(infiles), out_phy, (" (w/ flag)" if add_flag else ""))

        Node.__init__(self,
                      description  = description,
                      input_files  = infiles,
                      output_files = [out_phy],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        msa = MSA.join(*(MSA.from_file(filename) for filename in sorted(self.input_files)))

        with open(reroot_path(temp, self._out_phy), "w") as output:
            output.write(interleaved_phy(msa, add_flag = self._add_flag))


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_phy), self._out_phy)



class FastaToSequentialPhyNode(Node):
    def __init__(self, infiles, out_phy, add_flag = False, dependencies = ()):
        self._add_flag  = add_flag
        self._out_phy   = out_phy

        description  = "<FastaToInterleavedPhy: %i file(s) -> '%s'%s>" % \
            (len(infiles), out_phy, (" (w/ flag)" if add_flag else ""))

        Node.__init__(self,
                      description  = description,
                      input_files  = infiles,
                      output_files = [out_phy],
                      dependencies = dependencies)


    def _run(self, _config, temp):
        # Read and check that MSAs share groups
        msas = [MSA.from_file(filename) for filename in sorted(self.input_files)]
        MSA.validate(*msas)

        blocks = []
        for msa in msas:
            blocks.append(sequential_phy(msa, add_flag = self._add_flag))

        with open(reroot_path(temp, self._out_phy), "w") as output:
            output.write("\n\n".join(blocks))


    def  _teardown(self, _config, temp):
        move_file(reroot_path(temp, self._out_phy), self._out_phy)





def _read_sequences(filenames):
    results = {}
    for filename in filenames:
        results[filename] = MSA.from_file(filename)
    MSA.validate(*results.values())

    return results.iteritems()
