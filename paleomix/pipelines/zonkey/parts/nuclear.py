#!/usr/bin/python3
#
# Copyright (c) 2016 Mikkel Schubert <MikkelSch@gmail.com>
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
import hashlib
import itertools
import math
import os
import random
from typing import Iterable, Optional

import pysam

import paleomix.common.fileutils as fileutils
import paleomix.common.rtools as rtools
import paleomix.common.versions as versions
import paleomix.tools.factory as factory
from paleomix.common.command import (
    AtomicCmd,
    AuxiliaryFile,
    InputFile,
    OutputFile,
    SequentialCmds,
    TempOutputFile,
)
from paleomix.node import CommandNode, Node, NodeError
from paleomix.pipelines.zonkey.common import RSCRIPT_VERSION, read_summary

ADMIXTURE_VERSION = versions.Requirement(
    call=("admixture", "--version"),
    regexp=r"(\d+\.\d+)",
    specifiers=">=1.3",
)

PLINK_VERSION = versions.Requirement(
    call=("plink", "--noweb", "--help", "--out", "/tmp/plink"),
    regexp=r"v(\d+\.\d+)",
    specifiers=">=1.7",
)

SMARTPCA_VERSION = versions.Requirement(
    call=("smartpca",),
    regexp=r"version: (\d+)",
    specifiers=">=13050",
)

TREEMIX_VERSION = versions.Requirement(
    call=("treemix",),
    regexp=r"TreeMix v. (\d+\.\d+)",
    specifiers=">=1.12",
)


class BuildTPEDFilesNode(CommandNode):
    def __init__(self, output_root, table, bamfile, downsample, dependencies=()):
        command = factory.new(
            [
                "zonkey:tped",
                "--name",
                "Sample",
                "--downsample",
                downsample,
                "%(TEMP_DIR)s",
                InputFile(table),
                InputFile(bamfile),
            ],
            extra_files=[
                OutputFile(os.path.join(output_root, "common.tfam")),
                OutputFile(os.path.join(output_root, "common.summary")),
                OutputFile(os.path.join(output_root, "incl_ts.tped")),
                OutputFile(os.path.join(output_root, "excl_ts.tped")),
            ],
        )

        if not downsample:
            # Needed for random access (chromosomes are read 1 ... 31)
            command.add_extra_files([InputFile(bamfile + ".bai")])

        CommandNode.__init__(
            self,
            description="building TPED file from %s" % (bamfile,),
            command=command,
            dependencies=dependencies,
        )


class BuildBEDFilesNode(CommandNode):
    def __init__(
        self,
        output_prefix: str,
        tfam: str,
        tped: str,
        plink_parameters=None,
        dependencies: Iterable[Node] = (),
    ):
        temp_prefix = os.path.basename(output_prefix)

        command = AtomicCmd(
            [
                "plink",
                "--make-bed",
                "--noweb",
                "--tped",
                InputFile(tped),
                "--tfam",
                InputFile(tfam),
                "--out",
                TempOutputFile(temp_prefix),
            ]
            + self._parse_parameters(plink_parameters),
            extra_files=[
                OutputFile(output_prefix + ".bed"),
                OutputFile(output_prefix + ".bim"),
                OutputFile(output_prefix + ".fam"),
                OutputFile(output_prefix + ".log"),
                TempOutputFile(temp_prefix + ".nosex"),
                TempOutputFile(temp_prefix + ".nof"),
            ],
            requirements=[PLINK_VERSION],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="building BED files from %s" % (tped,),
            command=command,
            dependencies=dependencies,
        )

    @classmethod
    def _parse_parameters(cls, parameters):
        return parameters.split()


class AdmixtureNode(CommandNode):
    def __init__(
        self,
        input_file: str,
        k_groups: int,
        output_root: str,
        groups,
        dependencies: Iterable[Node] = (),
    ):
        self._groups = groups
        self._input_file = input_file

        prefix = os.path.splitext(os.path.basename(input_file))[0]
        output_prefix = os.path.join(output_root, "%s.%i" % (prefix, k_groups))

        command = AtomicCmd(
            [
                "admixture",
                "-s",
                random.randint(0, 2 ** 16 - 1),
                "--supervised",
                TempOutputFile(prefix + ".bed"),
                int(k_groups),
            ],
            stdout=output_prefix + ".log",
            extra_files=[
                InputFile(input_file),
                InputFile(fileutils.swap_ext(input_file, ".bim")),
                InputFile(fileutils.swap_ext(input_file, ".fam")),
                TempOutputFile(prefix + ".bim"),
                TempOutputFile(prefix + ".fam"),
                TempOutputFile(prefix + ".pop"),
                OutputFile(output_prefix + ".P"),
                OutputFile(output_prefix + ".Q"),
            ],
            requirements=[ADMIXTURE_VERSION],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="estimating admixture from %s.*" % (input_file,),
            command=command,
            dependencies=dependencies,
        )

    def _setup(self, temp):
        CommandNode._setup(self, temp)

        input_files = [
            self._input_file,
            fileutils.swap_ext(self._input_file, ".bim"),
            fileutils.swap_ext(self._input_file, ".fam"),
        ]

        for filename in input_files:
            basename = os.path.basename(filename)
            os.symlink(os.path.abspath(filename), os.path.join(temp, basename))

        fam_filename = fileutils.swap_ext(self._input_file, ".fam")
        pop_filename = fileutils.swap_ext(fam_filename, ".pop")
        pop_filename = fileutils.reroot_path(temp, pop_filename)

        with open(fam_filename) as fam_handle:
            with open(pop_filename, "w") as pop_handle:
                for line in fam_handle:
                    sample, _ = line.split(None, 1)

                    pop_handle.write("%s\n" % (self._groups.get(sample, "-"),))


class SelectBestAdmixtureNode(Node):
    def __init__(
        self,
        replicates: Iterable[AdmixtureNode],
        output_root: str,
        dependencies: Iterable[Node] = (),
    ):
        replicates = tuple(replicates)
        if not replicates:
            raise ValueError("No replicates passed to SelectBestAdmixture")

        input_files = []
        ref_filenames = None
        for node in replicates:
            filenames = frozenset(
                os.path.basename(filename) for filename in node.output_files
            )

            if ref_filenames is None:
                ref_filenames = filenames
            elif ref_filenames != filenames:
                raise RuntimeError(
                    "Node %r does not contain expected files, "
                    "%r, vs %r" % (node, ref_filenames, filenames)
                )

            input_files.extend(node.output_files)

        assert ref_filenames is not None
        output_files = [
            os.path.join(output_root, filename) for filename in ref_filenames
        ]

        self._files = tuple(node.output_files for node in replicates)
        self._output_root = output_root

        Node.__init__(
            self,
            description="selecting best admixture result from %s" % (output_root,),
            input_files=input_files,
            output_files=output_files,
            dependencies=tuple(dependencies) + tuple(replicates),
        )

    def _run(self, temp):
        likelihoods = []
        for fileset in self._files:
            for filename in fileset:
                if filename.endswith(".log"):
                    likelihoods.append((self._read_admixture_log(filename), fileset))
                    break
            else:
                raise NodeError(
                    "No log-file found in list of admixture "
                    "output-files: %r" % (fileset,)
                )

        _, fileset = max(likelihoods)
        for src_filename in fileset:
            dst_filename = fileutils.reroot_path(self._output_root, src_filename)
            fileutils.copy_file(src_filename, dst_filename)

    @classmethod
    def _read_admixture_log(cls, filename):
        with open(filename) as handle:
            for line in handle:
                if line.startswith("Loglikelihood:"):
                    return float(line.split()[1])

            raise NodeError(
                "Could not find likelihood value in log-file %r; "
                "looking for line starting with 'Loglikelihood:'" % (filename,)
            )


class AdmixturePlotNode(CommandNode):
    def __init__(self, input_file, output_prefix, order, samples, dependencies=()):
        self._samples = samples
        self._order = tuple(order) + ("Sample",)

        command = AtomicCmd(
            (
                "Rscript",
                AuxiliaryFile(rtools.rscript("zonkey", "admixture.r")),
                InputFile(input_file),
                TempOutputFile("samples.txt"),
                TempOutputFile(output_prefix),
            ),
            extra_files=[
                InputFile(samples),
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".png"),
            ],
            requirements=[
                RSCRIPT_VERSION,
                rtools.requirement("ggplot2"),
                rtools.requirement("reshape2"),
            ],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="plotting admixture results from %s" % (input_file,),
            command=command,
            dependencies=dependencies,
        )

    def _setup(self, temp):
        samples = {}
        with open(self._samples) as handle:
            header = handle.readline().strip().split("\t")
            for line in handle:
                row = dict(zip(header, line.strip().split("\t")))
                samples[row["Name"]] = row

        with open(os.path.join(temp, "samples.txt"), "w") as handle:
            handle.write("{}\n".format("\t".join(header)))

            for name in self._order:
                row = samples[name]
                handle.write("{}\n".format("\t".join(row[key] for key in header)))

        CommandNode._setup(self, temp)


class BuildFreqFilesNode(CommandNode):
    def __init__(
        self,
        input_prefix: str,
        output_prefix: str,
        tfam: str,
        parameters: Optional[str] = None,
        dependencies: Iterable[Node] = (),
    ):
        basename = os.path.basename(output_prefix)

        plink_cmd = [
            "plink",
            "--freq",
            "--missing",
            "--noweb",
            "--bfile",
            os.path.abspath(input_prefix),
            "--within",
            TempOutputFile("samples.clust"),
            "--out",
            TempOutputFile(basename),
        ]

        if parameters:
            plink_cmd.extend(parameters.split())

        plink = AtomicCmd(
            plink_cmd,
            extra_files=[
                InputFile(input_prefix + ".bed"),
                InputFile(input_prefix + ".bim"),
                InputFile(input_prefix + ".fam"),
                OutputFile(output_prefix + ".frq.strat.nosex"),
                OutputFile(output_prefix + ".frq.strat.log"),
                TempOutputFile(basename + ".imiss"),
                TempOutputFile(basename + ".lmiss"),
            ],
            requirements=[
                PLINK_VERSION,
            ],
            set_cwd=True,
        )

        gzip = AtomicCmd(
            ["gzip", TempOutputFile(basename + ".frq.strat")],
            extra_files=[
                OutputFile(output_prefix + ".frq.strat.gz"),
            ],
        )

        self._tfam = tfam
        self._basename = basename
        CommandNode.__init__(
            self,
            description="building frequency files from %s.*" % (input_prefix,),
            command=SequentialCmds((plink, gzip)),
            dependencies=dependencies,
        )

    def _setup(self, temp):
        CommandNode._setup(self, temp)

        with open(self._tfam) as in_handle:
            samples = [line.split(None, 1)[0] for line in in_handle]

        with open(os.path.join(temp, "samples.clust"), "w") as handle:
            for sample in samples:
                handle.write("{0} {0} {0}\n".format(sample))

    def _teardown(self, temp):
        for ext in ("log", "nosex"):
            os.rename(
                os.path.join(temp, self._basename + "." + ext),
                os.path.join(temp, self._basename + ".frq.strat." + ext),
            )

        CommandNode._teardown(self, temp)


class FreqToTreemixNode(Node):
    def __init__(self, input_file, output_file, dependencies=()):
        Node.__init__(
            self,
            description="converting %s to treemix format" % (input_file,),
            input_files=(input_file,),
            output_files=(output_file,),
            dependencies=dependencies,
        )

    def _run(self, temp):
        (input_file,) = self.input_files
        (output_file,) = self.output_files
        temp_filename = os.path.basename(output_file)

        with open(os.path.join(temp, temp_filename), "w") as handle:
            header = None
            table = self._parse_freq_table(input_file)

            for _, rows in itertools.groupby(table, lambda row: row[:2]):
                if header is None:
                    rows = tuple(rows)  # Must not consume iterator
                    header = list(sorted(row[2] for row in rows))
                    handle.write("%s\n" % (" ".join(header)))

                result = []
                rows = dict((row[2], row) for row in rows)
                for sample in header:
                    _, _, _, mac, nchroms = rows[sample]
                    result.append("%s,%i" % (mac, int(nchroms) - int(mac)))

                handle.write("%s\n" % (" ".join(result),))

    def _teardown(self, temp):
        (output_file,) = self.output_files
        temp_file = os.path.join(temp, os.path.basename(output_file))

        fileutils.move_file(temp_file, output_file)
        Node._teardown(self, temp)

    @classmethod
    def _parse_freq_table(cls, filename):
        with fileutils.open_rt(filename) as handle:
            handle.readline()  # Skip header

            for line in handle:
                chrom, snp, clst, _, _, _, mac, nchroms = line.split()

                yield (chrom, snp, clst, int(mac), int(nchroms))


class TreemixNode(CommandNode):
    def __init__(
        self,
        data,
        input_file: str,
        output_prefix: str,
        m: int = 0,
        k=100,
        outgroup=(),
        dependencies: Iterable[Node] = (),
    ):
        self._param_m = m
        self._param_outgroup = outgroup
        self._params_file = output_prefix + ".parameters.txt"
        self._parameters_hash = "%s.%s" % (
            output_prefix,
            hash_params(k=k, m=m, global_set=True, outgroup=tuple(sorted(outgroup))),
        )

        command = AtomicCmd(
            [
                "treemix",
                "-i",
                InputFile(input_file),
                "-o",
                TempOutputFile(output_prefix),
                "-global",
                "-m",
                m,
            ],
            extra_files=[
                OutputFile(output_prefix + ".cov.gz"),
                OutputFile(output_prefix + ".covse.gz"),
                OutputFile(output_prefix + ".edges.gz"),
                OutputFile(output_prefix + ".llik"),
                OutputFile(output_prefix + ".modelcov.gz"),
                OutputFile(output_prefix + ".treeout.gz"),
                OutputFile(output_prefix + ".vertices.gz"),
                OutputFile(self._params_file),
                OutputFile(self._parameters_hash),
            ],
            requirements=[TREEMIX_VERSION],
            set_cwd=True,
        )

        if outgroup:
            command.append("-root", ",".join(outgroup))

        if isinstance(k, int):
            command.append("-k", k)
            self._param_k = k
            self._k_file = self._k_field = None
        elif isinstance(k, tuple) and all(isinstance(v, str) for v in k):
            self._k_field, self._k_file = k
            self._genome_size = sum(value["Size"] for value in data.contigs.values())
            self._snp_distance = data.settings["SNPDistance"]
        else:
            raise ValueError("k must be int or (key, path) in TreemixNode")

        CommandNode.__init__(
            self,
            description="running treemix on %s.*" % (input_file,),
            command=command,
            dependencies=dependencies,
        )

    def _setup(self, temp: str):
        if self._k_file is not None:
            stats = read_summary(self._k_file)
            n_sites = float(stats[self._k_field])
            k = max(
                1, int(math.ceil(self._snp_distance / (self._genome_size / n_sites)))
            )

            self._param_k = k
            assert isinstance(self._command, AtomicCmd)
            self._command._command.extend(("-k", str(k)))

        CommandNode._setup(self, temp)

    def _teardown(self, temp: str):
        with open(fileutils.reroot_path(temp, self._params_file), "w") as out:
            out.write("k: %i\n" % (self._param_k,))
            out.write("m: %i\n" % (self._param_m,))
            out.write("outgroup: %r\n" % (list(self._param_outgroup),))

        open(fileutils.reroot_path(temp, self._parameters_hash), "w").close()

        CommandNode._teardown(self, temp)


class PlotTreemixNode(CommandNode):
    def __init__(
        self,
        samples: str,
        prefix: str,
        output_prefix: str,
        dependencies: Iterable[Node] = (),
    ):
        abs_prefix = os.path.abspath(prefix)
        basename = os.path.basename(output_prefix)

        # TreeMix plots with migration edges
        cmd_1 = self._plot_command(
            input_prefix=prefix,
            args=(
                "plot_tree",
                abs_prefix,
                InputFile(samples),
                TempOutputFile(basename + "_tree"),
            ),
            extra_files=[
                OutputFile(output_prefix + "_tree.pdf"),
                OutputFile(output_prefix + "_tree.png"),
            ],
        )

        # Heatmap showing TreeMix residuals
        cmd_2 = self._plot_command(
            input_prefix=prefix,
            args=(
                "plot_residuals",
                abs_prefix,
                InputFile(samples),
                TempOutputFile(basename + "_residuals"),
            ),
            extra_files=[
                OutputFile(output_prefix + "_residuals.pdf"),
                OutputFile(output_prefix + "_residuals.png"),
            ],
        )

        # Text file containing % of variance explained by model
        cmd_3 = self._plot_command(
            input_prefix=prefix,
            args=(
                "variance",
                abs_prefix,
                OutputFile(output_prefix + "_variance.txt"),
            ),
        )

        CommandNode.__init__(
            self,
            description="plotting treemix results to %s.*" % output_prefix,
            command=SequentialCmds((cmd_1, cmd_2, cmd_3)),
            dependencies=dependencies,
        )

    @classmethod
    def _plot_command(cls, input_prefix, args, extra_files=[], **kwargs):
        script = rtools.rscript("zonkey", "treemix.r")

        return AtomicCmd(
            ("Rscript", AuxiliaryFile(script)) + tuple(args),
            extra_files=[
                InputFile(input_prefix + ".cov.gz"),
                InputFile(input_prefix + ".covse.gz"),
                InputFile(input_prefix + ".edges.gz"),
                InputFile(input_prefix + ".modelcov.gz"),
                InputFile(input_prefix + ".vertices.gz"),
            ]
            + extra_files,
            requirements=[
                RSCRIPT_VERSION,
                rtools.requirement("RColorBrewer"),
            ],
            set_cwd=True,
            **kwargs,
        )


class SmartPCANode(CommandNode):
    def __init__(
        self,
        input_prefix: str,
        output_prefix: str,
        nchroms: int,
        dependencies: Iterable[Node] = (),
    ):
        self._input_prefix = input_prefix
        self._output_prefix = output_prefix
        self._nchroms = nchroms

        cmd = AtomicCmd(
            ("smartpca", "-p", TempOutputFile("parameters.txt")),
            stdout=output_prefix + ".log",
            extra_files=[
                InputFile(input_prefix + ".bed"),
                InputFile(input_prefix + ".bim"),
                InputFile(input_prefix + ".fam"),
                OutputFile(output_prefix + ".evec"),
                OutputFile(output_prefix + ".eval"),
                OutputFile(output_prefix + ".deleted_snps"),
            ],
            requirements=[SMARTPCA_VERSION],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="performing PCA on %s" % (input_prefix,),
            command=cmd,
            dependencies=dependencies,
        )

    def _setup(self, temp):
        CommandNode._setup(self, temp)

        with open(os.path.join(temp, "parameters.txt"), "w") as handle:
            handle.write(
                """
genotypename:      {input_prefix}.bed
snpname:           {input_prefix}.bim
indivname:         {input_prefix}.fam
evecoutname:       {output_prefix}.evec
evaloutname:       {output_prefix}.eval
deletsnpoutname:   {output_prefix}.deleted_snps
altnormstyle:      NO
numoutevec:        5
familynames:       YES
numoutlieriter:    1
numchrom:          {nchroms}
numthreads:        1
""".format(
                    input_prefix=os.path.abspath(self._input_prefix),
                    output_prefix=os.path.basename(self._output_prefix),
                    nchroms=self._nchroms,
                )
            )

    def _teardown(self, temp):
        # Ensure that this file exists even when no filtered SNPs
        deleted_snps = os.path.basename(self._output_prefix) + ".deleted_snps"
        open(os.path.join(temp, deleted_snps), "a").close()

        CommandNode._teardown(self, temp)


class PlotPCANode(CommandNode):
    def __init__(self, samples, prefix, output_prefix, dependencies=()):
        cmd = AtomicCmd(
            [
                "Rscript",
                AuxiliaryFile(rtools.rscript("zonkey", "pca.r")),
                os.path.abspath(prefix),
                InputFile(samples),
                TempOutputFile(output_prefix),
            ],
            extra_files=[
                InputFile(prefix + ".eval"),
                InputFile(prefix + ".evec"),
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".png"),
            ],
            requirements=[
                RSCRIPT_VERSION,
                rtools.requirement("ggplot2"),
                rtools.requirement("ggrepel"),
            ],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="plotting PCA to %s.*" % (output_prefix,),
            command=cmd,
            dependencies=dependencies,
        )


class PlotCoverageNode(CommandNode):
    def __init__(self, contigs, mapping, input_file, output_prefix, dependencies=()):
        self._contigs = contigs
        self._mapping = dict(zip(mapping.values(), mapping))
        self._input_file = input_file

        cmd = AtomicCmd(
            (
                "Rscript",
                AuxiliaryFile(rtools.rscript("zonkey", "coverage.r")),
                TempOutputFile("contigs.table"),
                TempOutputFile(output_prefix),
            ),
            extra_files=[
                InputFile(input_file),
                OutputFile(output_prefix + ".pdf"),
                OutputFile(output_prefix + ".png"),
            ],
            requirements=[
                RSCRIPT_VERSION,
                rtools.requirement("ggplot2"),
            ],
            set_cwd=True,
        )

        CommandNode.__init__(
            self,
            description="plotting coverage to %s.*" % (output_prefix,),
            command=cmd,
            dependencies=dependencies,
        )

    def _setup(self, temp):
        with open(os.path.join(temp, "contigs.table"), "w") as handle:
            handle.write("ID\tSize\tNs\tHits\n")

            # Workaround for pysam < 0.9 returning list, >= 0.9 returning str
            for line in "".join(pysam.idxstats(self._input_file)).split("\n"):
                line = line.strip()
                if not line:
                    continue

                name, size, hits, _ = line.split("\t")
                name = self._mapping.get(name, name)
                if name not in self._contigs:
                    # Excluding contigs is allowed
                    continue

                row = {
                    "ID": name,
                    "Size": self._contigs[name]["Size"],
                    "Ns": self._contigs[name]["Ns"],
                    "Hits": hits,
                }

                handle.write("{ID}\t{Size}\t{Ns}\t{Hits}\n".format(**row))

        CommandNode._setup(self, temp)


def hash_params(*args, **kwargs):
    return hashlib.md5(repr([args, kwargs]).encode("utf-8")).hexdigest()
