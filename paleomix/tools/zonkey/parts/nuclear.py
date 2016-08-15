#!/usr/bin/python
#
# Copyright (c) 2016 Mikkel Schubert <MSchubert@snm.ku.dk>
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

import pysam

import paleomix.common.fileutils as fileutils
import paleomix.common.rtools as rtools
import paleomix.common.versions as versions
import paleomix.tools.factory as factory

from paleomix.atomiccmd.builder import AtomicCmdBuilder
from paleomix.atomiccmd.command import AtomicCmd
from paleomix.atomiccmd.sets import SequentialCmds
from paleomix.node import CommandNode, Node, NodeError

from paleomix.tools.zonkey.common import \
    RSCRIPT_VERSION, \
    contig_name_to_plink_name, \
    read_summary


ADMIXTURE_VERSION = versions.Requirement(call=("admixture", "--version"),
                                         search="(\d+)\.(\d+)",
                                         checks=versions.GE(1, 3))

PLINK_VERSION = versions.Requirement(call=("plink", "--noweb", "--help",
                                           "--out", "/tmp/plink"),
                                     search="v(\d+)\.(\d+)",
                                     checks=versions.GE(1, 7))

SMARTPCA_VERSION = versions.Requirement(call=("smartpca",),
                                        search="version: (\d+)",
                                        checks=versions.GE(13050))

TREEMIX_VERSION = versions.Requirement(call=("treemix",),
                                       search="TreeMix v. (\d+)\.(\d+)",
                                       checks=versions.GE(1, 12))


class BuildTPEDFilesNode(CommandNode):
    def __init__(self, output_root, table, bamfile, downsample,
                 dependencies=()):
        cmd = factory.new("zonkey_tped")
        cmd.set_option("--name", "Sample")
        cmd.set_option("--downsample", downsample)
        cmd.add_value("%(TEMP_DIR)s")
        cmd.add_value("%(IN_TABLE)s")
        cmd.add_value("%(IN_BAM)s")

        if not downsample:
            # Needed for random access (chromosomes are read 1 ... 31)
            cmd.set_kwargs(IN_BAI=fileutils.swap_ext(bamfile, ".bai"))

        cmd.set_kwargs(OUT_TFAM=os.path.join(output_root, "common.tfam"),
                       OUT_SUMMARY=os.path.join(output_root, "common.summary"),
                       OUT_TPED_INCL_TS=os.path.join(output_root,
                                                     "incl_ts.tped"),
                       OUT_TPED_EXCL_TS=os.path.join(output_root,
                                                     "excl_ts.tped"),
                       IN_TABLE=table,
                       IN_BAM=bamfile)

        CommandNode.__init__(self,
                             description="<BuildTPEDFiles -> %r>"
                             % (os.path.join(output_root, '*'),),
                             command=cmd.finalize(),
                             dependencies=dependencies)


class BuildBEDFilesNode(CommandNode):
    def __init__(self, output_prefix, tfam, tped,
                 indep_filter=None, indep_parameters=None,
                 plink_parameters=None,
                 dependencies=()):
        temp_prefix = os.path.basename(output_prefix)

        plink_cmd = ["plink", "--make-bed", "--noweb",
                     "--tped", "%(IN_TPED)s",
                     "--tfam", "%(IN_TFAM)s",
                     "--out", "%(TEMP_OUT_PREFIX)s"]

        plink_cmd.extend(self._parse_parameters(plink_parameters))

        command = AtomicCmd(plink_cmd,
                            IN_TPED=tped,
                            IN_TFAM=tfam,
                            TEMP_OUT_PREFIX=temp_prefix,
                            OUT_BED=output_prefix + ".bed",
                            OUT_BIM=output_prefix + ".bim",
                            OUT_FAM=output_prefix + ".fam",
                            OUT_LOG=output_prefix + ".log",
                            TEMP_OUT_NOSEX=temp_prefix + ".nosex",
                            TEMP_OUT_NOF=temp_prefix + ".nof",
                            CHECK_VERSION=PLINK_VERSION,
                            set_cwd=True)

        CommandNode.__init__(self,
                             description="<BuildBEDFiles -> '%s.*'>"
                             % (output_prefix,),
                             command=command,
                             dependencies=dependencies)

    @classmethod
    def _parse_parameters(cls, parameters):
        return parameters.split()


class BuildFilteredBEDFilesNode(CommandNode):
    def __init__(self, output_prefix, tfam, tped,
                 indep_filter=None, indep_parameters=None,
                 plink_parameters=None,
                 dependencies=()):

        assert indep_filter in ('indep',
                                'indep-pairphase',
                                'indep-pairwise'), indep_filter
        assert len(indep_parameters) == 3, indep_parameters

        parameters = self._parse_parameters(plink_parameters)

        plink_cmd = ["plink", "--noweb",
                     "--tped", "%(IN_TPED)s",
                     "--tfam", "%(IN_TFAM)s",
                     "--out", "%(TEMP_OUT_PREFIX)s",
                     '--' + indep_filter]
        plink_cmd.extend(indep_parameters)
        plink_cmd.extend(parameters)

        cmd_indep = AtomicCmd(plink_cmd,
                              IN_TFAM=tfam,
                              IN_TPED=tped,
                              TEMP_OUT_PREFIX="indep",
                              TEMP_OUT_LOG="indep.log",
                              TEMP_OUT_NOSEX="indep.nosex",
                              TEMP_OUT_PRUNE_IN="indep.prune.in",
                              TEMP_OUT_PRUNE_OUT="indep.prune.out",
                              set_cwd=True)

        basename = os.path.basename(output_prefix)
        cmd_filter = AtomicCmd(["plink", "--noweb", "--make-bed",
                                "--tped", "%(IN_TPED)s",
                                "--tfam", "%(IN_TFAM)s",
                                "--extract", "%(TEMP_IN_PRUNE)s",
                                "--out", "%(TEMP_OUT_PREFIX)s"] +
                               parameters,
                               IN_TFAM=tfam,
                               IN_TPED=tped,
                               TEMP_OUT_PREFIX=basename,
                               TEMP_IN_PRUNE="indep.prune.in",
                               TEMP_OUT_NOSEX=basename + ".nosex",
                               TEMP_OUT_LOG=basename + ".log",
                               OUT_LOG=output_prefix + ".log",
                               OUT_BED=output_prefix + ".bed",
                               OUT_BIM=output_prefix + ".bim",
                               OUT_FAM=output_prefix + ".fam",
                               set_cwd=True)

        CommandNode.__init__(self,
                             description="<BuildFilteredBEDFiles -> '%s.*'>"
                             % (output_prefix,),
                             command=SequentialCmds((cmd_indep, cmd_filter)),
                             dependencies=dependencies)

    @classmethod
    def _parse_parameters(cls, parameters):
        return parameters.split()


class AdmixtureNode(CommandNode):
    def __init__(self, input_file, k_groups, output_root,
                 samples=None, dependencies=()):
        self._samples = samples
        self._input_file = input_file
        self._k_groups = k_groups

        group_key = "Group(%i)" % (self._k_groups,)
        self._supervised = samples and any((row[group_key] != '-')
                                           for row in samples.itervalues())

        assert k_groups in (2, 3), k_groups
        prefix = os.path.splitext(os.path.basename(input_file))[0]
        output_prefix = os.path.join(output_root,
                                     "%s.%i" % (prefix, k_groups))

        cmd = AtomicCmdBuilder("admixture",
                               IN_FILE_BED=input_file,
                               IN_FILE_BIM=fileutils.swap_ext(input_file,
                                                              ".bim"),
                               IN_FILE_FAM=fileutils.swap_ext(input_file,
                                                              ".fam"),

                               TEMP_OUT_FILE_BED=prefix + ".bed",
                               TEMP_OUT_FILE_BIM=prefix + ".bim",
                               TEMP_OUT_FILE_FAM=prefix + ".fam",
                               TEMP_OUT_FILE_POP=prefix + ".pop",

                               OUT_P=output_prefix + ".P",
                               OUT_Q=output_prefix + ".Q",
                               OUT_STDOUT=output_prefix + ".log",

                               CHECK_VERSION=ADMIXTURE_VERSION,
                               set_cwd=True)

        cmd.set_option("-s", random.randint(0, 2 ** 16 - 1))

        if self._supervised:
            cmd.set_option("--supervised")

        cmd.add_value("%(TEMP_OUT_FILE_BED)s")
        cmd.add_value(int(k_groups))

        CommandNode.__init__(self,
                             description="<Admixture -> '%s.*''>"
                             % (output_prefix,),
                             command=cmd.finalize(),
                             dependencies=dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        input_files = [
            self._input_file,
            fileutils.swap_ext(self._input_file, ".bim"),
            fileutils.swap_ext(self._input_file, ".fam"),
        ]

        for filename in input_files:
            basename = os.path.basename(filename)
            os.symlink(os.path.abspath(filename), os.path.join(temp, basename))

        if self._supervised:
            fam_filename = fileutils.swap_ext(self._input_file, ".fam")

            pop_filename = fileutils.swap_ext(fam_filename, ".pop")
            pop_filename = fileutils.reroot_path(temp, pop_filename)

            key = "Group(%i)" % (self._k_groups,)
            with open(fam_filename) as fam_handle:
                with open(pop_filename, "w") as pop_handle:
                    for line in fam_handle:
                        sample, _ = line.split(None, 1)
                        group = self._samples.get(sample, {}).get(key, "-")

                        pop_handle.write("%s\n" % (group,))


class SelectBestAdmixtureNode(Node):
    def __init__(self, replicates, output_root, dependencies=()):
        replicates = tuple(replicates)
        if not replicates:
            raise ValueError("No replicates passed to SelectBestAdmixture")

        input_files = []
        ref_filenames = None
        for node in replicates:
            filenames = frozenset(os.path.basename(filename)
                                  for filename in node.output_files)

            if ref_filenames is None:
                ref_filenames = filenames
            elif ref_filenames != filenames:
                raise RuntimeError("Node %r does not contain expected files, "
                                   "%r, vs %r" % (node, ref_filenames,
                                                  filenames))

            input_files.extend(node.output_files)

        output_files = [os.path.join(output_root, filename)
                        for filename in ref_filenames]

        self._ref_filenames = ref_filenames
        self._files = tuple(node.output_files for node in replicates)
        self._output_root = output_root

        Node.__init__(self,
                      description="<SelectBestAdmixture -> %r>"
                      % (output_root,),
                      input_files=input_files,
                      output_files=output_files,
                      dependencies=tuple(dependencies) + tuple(replicates))

    def _run(self, config, temp):
        likelihoods = []
        for fileset in self._files:
            for filename in fileset:
                if filename.endswith(".log"):
                    likelihoods.append((self._read_admixture_log(filename),
                                        fileset))
                    break
            else:
                raise NodeError("No log-file found in list of admixture "
                                "output-files: %r" % (fileset,))

        _, fileset = max(likelihoods)
        for src_filename in fileset:
            dst_filename = fileutils.reroot_path(self._output_root,
                                                 src_filename)
            fileutils.copy_file(src_filename, dst_filename)

    @classmethod
    def _read_admixture_log(cls, filename):
        with open(filename) as handle:
            for line in handle:
                if line.startswith("Loglikelihood:"):
                    return float(line.split()[1])

            raise NodeError("Could not find likelihood value in log-file %r; "
                            "looking for line starting with 'Loglikelihood:'"
                            % (filename,))


class AdmixturePlotNode(CommandNode):
    def __init__(self, input_file, output_prefix, order, samples,
                 dependencies=()):
        self._samples = samples
        self._order = tuple(order) + ("Sample",)

        script = rtools.rscript("zonkey", "admixture.r")

        cmd = AtomicCmd(("Rscript", script, "%(IN_FILE)s",
                         "%(TEMP_OUT_NAMES)s", "%(TEMP_OUT_PREFIX)s"),
                        AUX_RSCRIPT=script,
                        IN_FILE=input_file,
                        IN_SAMPLES=samples,
                        OUT_PDF=output_prefix + ".pdf",
                        OUT_PNG=output_prefix + ".png",
                        TEMP_OUT_NAMES="samples.txt",
                        TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                        CHECK_R=RSCRIPT_VERSION,
                        CHECK_R_GGPLOT2=rtools.requirement("ggplot2"),
                        CHECK_R_RESHAPE2=rtools.requirement("reshape2"),
                        set_cwd=True)

        CommandNode.__init__(self,
                             description="<AdmixturePlot -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        samples = {}
        with open(self._samples) as handle:
            header = handle.readline().strip().split('\t')
            for line in handle:
                row = dict(zip(header, line.strip().split('\t')))
                samples[row["Name"]] = row

        with open(os.path.join(temp, "samples.txt"), "w") as handle:
            handle.write("{}\n".format("\t".join(header)))

            for name in self._order:
                row = samples[name]
                handle.write("{}\n".format("\t".join(row[key]
                                                     for key in header)))

        CommandNode._setup(self, config, temp)


class BuildFreqFilesNode(CommandNode):
    def __init__(self, input_prefix, output_prefix, tfam,
                 parameters=None, dependencies=()):
        basename = os.path.basename(output_prefix)

        plink_cmd = ["plink", "--freq", "--missing", "--noweb",
                     "--bfile", os.path.abspath(input_prefix),
                     "--within", "%(TEMP_OUT_CLUST)s",
                     "--out", "%(TEMP_OUT_PREFIX)s"]

        if parameters:
            plink_cmd.extend(parameters.split())

        plink = AtomicCmd(plink_cmd,
                          IN_BED=input_prefix + ".bed",
                          IN_BIM=input_prefix + ".bim",
                          IN_FAM=input_prefix + ".fam",
                          TEMP_OUT_CLUST="samples.clust",
                          OUT_NOSEX=output_prefix + ".frq.strat.nosex",
                          OUT_LOG=output_prefix + ".frq.strat.log",
                          TEMP_OUT_PREFIX=basename,
                          CHECK_VERSION=PLINK_VERSION,
                          set_cwd=True)

        gzip = AtomicCmd(["gzip", "%(TEMP_IN_FREQ)s"],
                         TEMP_IN_FREQ=basename + ".frq.strat",
                         OUT_FREQ=output_prefix + ".frq.strat.gz")

        self._tfam = tfam
        self._basename = basename
        CommandNode.__init__(self,
                             description="<BuildFreqFiles -> '%s.*'"
                             % (output_prefix,),
                             command=SequentialCmds((plink, gzip)),
                             dependencies=dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        with open(self._tfam) as in_handle:
            samples = [line.split(None, 1)[0] for line in in_handle]

        with open(os.path.join(temp, "samples.clust"), "w") as handle:
            for sample in samples:
                handle.write("{0} {0} {0}\n".format(sample))

    def _teardown(self, config, temp):
        for ext in ("log", "nosex"):
            os.rename(os.path.join(temp, self._basename + "." + ext),
                      os.path.join(temp, self._basename + ".frq.strat." + ext))

        CommandNode._teardown(self, config, temp)


class FreqToTreemixNode(Node):
    def __init__(self, input_file, output_file, dependencies=()):
        Node.__init__(self,
                      description="<FreqToTreemix -> %r" % (output_file,),
                      input_files=(input_file,),
                      output_files=(output_file,),
                      dependencies=dependencies)

    def _run(self, _config, temp):
        input_file, = self.input_files
        output_file, = self.output_files
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

    def _teardown(self, config, temp):
        output_file, = self.output_files
        temp_file = os.path.join(temp, os.path.basename(output_file))

        fileutils.move_file(temp_file, output_file)
        Node._teardown(self, config, temp)

    @classmethod
    def _parse_freq_table(cls, filename):
        with fileutils.open_ro(filename) as handle:
            handle.readline()  # Skip header

            for line in handle:
                chrom, snp, clst, _, _, _, mac, nchroms = line.split()

                yield (chrom, snp, clst, int(mac), int(nchroms))


class TreemixNode(CommandNode):
    def __init__(self, data, input_file, output_prefix, m=0, k=100,
                 outgroup=(), dependencies=()):
        call = ["treemix",
                "-i", "%(IN_FILE)s",
                "-o", "%(TEMP_OUT_PREFIX)s",
                "-global",
                "-m", m]

        if outgroup:
            call.extend(("-root", ",".join(outgroup)))

        self._param_m = m
        self._param_outgroup = outgroup
        self._params_file = output_prefix + ".parameters.txt"

        if isinstance(k, int):
            call.extend(("-k", k))
            self._param_k = k
            self._k_file = self._k_field = None
        elif isinstance(k, tuple) and all(isinstance(v, str) for v in k):
            self._k_field, self._k_file = k
            self._genome_size = sum(value["Size"]
                                    for value in data.contigs.itervalues())
            self._snp_distance = data.settings["SNPDistance"]
        else:
            raise ValueError("k must be int or (key, path) in TreemixNode")

        self._parameters_hash \
            = "%s.%s" % (output_prefix,
                         hash_params(k=k, m=m, global_set=True,
                                     outgroup=tuple(sorted(outgroup))))

        cmd = AtomicCmd(call,
                        IN_FILE=input_file,
                        TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                        OUT_FILE_COV=output_prefix + ".cov.gz",
                        OUT_FILE_COVSE=output_prefix + ".covse.gz",
                        OUT_FILE_EDGES=output_prefix + ".edges.gz",
                        OUT_FILE_LLIK=output_prefix + ".llik",
                        OUT_FILE_MODELCOV=output_prefix + ".modelcov.gz",
                        OUT_FILE_TREEOUT=output_prefix + ".treeout.gz",
                        OUT_FILE_VERTICES=output_prefix + ".vertices.gz",
                        OUT_FILE_PARAMS=self._params_file,
                        OUT_FILE_PARAMS_HASH=self._parameters_hash,
                        CHECK_VERSION=TREEMIX_VERSION,
                        set_cwd=True)

        CommandNode.__init__(self,
                             description="<Treemix -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        if self._k_file is not None:
            stats = read_summary(self._k_file)
            n_sites = float(stats[self._k_field])
            k = max(1, int(math.ceil(self._snp_distance /
                                     (self._genome_size / n_sites))))

            self._param_k = k
            self._command._command.extend(("-k", str(k)))

        CommandNode._setup(self, config, temp)

    def _teardown(self, config, temp):
        with open(fileutils.reroot_path(temp, self._params_file), "w") as out:
            out.write("k: %i\n" % (self._param_k,))
            out.write("m: %i\n" % (self._param_m,))
            out.write("outgroup: %r\n" % (list(self._param_outgroup),))

        open(fileutils.reroot_path(temp, self._parameters_hash), "w").close()

        CommandNode._teardown(self, config, temp)


class PlotTreemixNode(CommandNode):
    def __init__(self, samples, prefix, output_prefix, dependencies=()):
        abs_prefix = os.path.abspath(prefix)
        basename = os.path.basename(output_prefix)

        # TreeMix plots with migration edges
        cmd_1 = self._plot_command(prefix, "plot_tree", abs_prefix,
                                   "%(IN_SAMPLES)s", "%(TEMP_OUT_PREFIX)s",
                                   IN_SAMPLES=samples,
                                   TEMP_OUT_PREFIX=basename + "_tree",
                                   OUT_PDF=output_prefix + "_tree.pdf",
                                   OUT_PNG=output_prefix + "_tree.png")

        # Heatmap showing TreeMix residuals
        cmd_2 = self._plot_command(prefix, "plot_residuals", abs_prefix,
                                   "%(IN_SAMPLES)s", "%(TEMP_OUT_PREFIX)s",
                                   IN_SAMPLES=samples,
                                   TEMP_OUT_PREFIX=basename + "_residuals",
                                   OUT_PDF=output_prefix + "_residuals.pdf",
                                   OUT_PNG=output_prefix + "_residuals.png")

        # Text file containing % of variance explained by model
        cmd_3 = self._plot_command(prefix, "variance", abs_prefix,
                                   "%(OUT_TXT)s",
                                   OUT_TXT=output_prefix + "_variance.txt")

        CommandNode.__init__(self,
                             description="<PlotTreemix -> '%s.*'>"
                             % (output_prefix,),
                             command=SequentialCmds((cmd_1, cmd_2, cmd_3)),
                             dependencies=dependencies)

    @classmethod
    def _plot_command(cls, input_prefix, *args, **kwargs):
        script = rtools.rscript("zonkey", "treemix.r")

        return AtomicCmd(("Rscript", script) + args,
                         AUX_RSCRIPT=script,
                         IN_FILE_COV=input_prefix + ".cov.gz",
                         IN_FILE_COVSE=input_prefix + ".covse.gz",
                         IN_FILE_EDGES=input_prefix + ".edges.gz",
                         IN_FILE_MODELCOV=input_prefix + ".modelcov.gz",
                         IN_FILE_VERTICES=input_prefix + ".vertices.gz",
                         CHECK_R=RSCRIPT_VERSION,
                         CHECK_R_BREW=rtools.requirement("RColorBrewer"),
                         set_cwd=True,
                         **kwargs)


class SmartPCANode(CommandNode):
    def __init__(self, input_prefix, output_prefix, nchroms, dependencies=()):
        self._input_prefix = input_prefix
        self._output_prefix = output_prefix
        self._nchroms = nchroms

        cmd = AtomicCmd(("smartpca", "-p", "%(TEMP_OUT_PARAMS)s"),
                        TEMP_OUT_PARAMS="parameters.txt",
                        IN_FILE_BED=input_prefix + ".bed",
                        IN_FILE_BIM=input_prefix + ".bim",
                        IN_FILE_FAM=input_prefix + ".fam",
                        OUT_STDOUT=output_prefix + ".log",
                        OUT_EVEC=output_prefix + ".evec",
                        OUT_EVAL=output_prefix + ".eval",
                        OUT_SNPS=output_prefix + ".deleted_snps",
                        CHECK_VERSION=SMARTPCA_VERSION,
                        set_cwd=True)

        CommandNode.__init__(self,
                             description="<SmartPCA -> '%s.*>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        CommandNode._setup(self, config, temp)

        with open(os.path.join(temp, "parameters.txt"), "w") as handle:
            handle.write("""
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
""".format(input_prefix=os.path.abspath(self._input_prefix),
           output_prefix=os.path.basename(self._output_prefix),
           nchroms=self._nchroms))

    def _teardown(self, config, temp):
        # Ensure that this file exists even when no filtered SNPs
        deleted_snps = os.path.basename(self._output_prefix) + ".deleted_snps"
        open(os.path.join(temp, deleted_snps), "a").close()

        CommandNode._teardown(self, config, temp)


class PlotPCANode(CommandNode):
    def __init__(self, samples, prefix, output_prefix, dependencies=()):
        abs_prefix = os.path.abspath(prefix)

        script = rtools.rscript("zonkey", "pca.r")
        call = ["Rscript", script,
                abs_prefix, "%(IN_SAMPLES)s", "%(TEMP_OUT_PREFIX)s"]

        cmd = AtomicCmd(call,
                        AUX_SCRIPT=script,
                        IN_FILE_EVAL=prefix + ".eval",
                        IN_FILE_EVEC=prefix + ".evec",
                        IN_SAMPLES=samples,
                        TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                        OUT_PDF=output_prefix + ".pdf",
                        OUT_PNG=output_prefix + ".png",
                        CHECK_R=RSCRIPT_VERSION,
                        CHECK_R_GGPLOT2=rtools.requirement("ggplot2"),
                        CHECK_R_LABELS=rtools.requirement("ggrepel"),
                        set_cwd=True)

        CommandNode.__init__(self,
                             description="<PlotPCA -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)


class PlotCoverageNode(CommandNode):
    def __init__(self, contigs, input_file, output_prefix, dependencies=()):
        self._contigs = contigs
        self._input_file = input_file

        script = rtools.rscript("zonkey", "coverage.r")
        cmd = AtomicCmd(("Rscript", script,
                         "%(TEMP_OUT_TABLE)s", "%(TEMP_OUT_PREFIX)s"),
                        AUX_RSCRIPT=script,
                        IN_FILE=input_file,
                        TEMP_OUT_TABLE="contigs.table",
                        OUT_PDF=output_prefix + ".pdf",
                        OUT_PNG=output_prefix + ".png",
                        TEMP_OUT_PREFIX=os.path.basename(output_prefix),
                        CHECK_R=RSCRIPT_VERSION,
                        CHECK_R_GGPLOT2=rtools.requirement("ggplot2"),
                        set_cwd=True)

        CommandNode.__init__(self,
                             description="<CoveragePlot -> '%s.*'>"
                             % (output_prefix,),
                             command=cmd,
                             dependencies=dependencies)

    def _setup(self, config, temp):
        with open(os.path.join(temp, "contigs.table"), "w") as handle:
            handle.write("ID\tSize\tNs\tHits\n")

            # Workaround for pysam < 0.9 returning list, >= 0.9 returning str
            for line in "".join(pysam.idxstats(self._input_file)).split('\n'):
                line = line.strip()
                if not line:
                    continue

                name, size, hits, _ = line.split('\t')
                name = contig_name_to_plink_name(name)
                if name is None or not (name.isdigit() or name == 'X'):
                    continue
                elif name not in self._contigs:
                    # Excluding contigs is allowed
                    continue

                if int(size) != self._contigs[name]['Size']:
                    raise NodeError("Size mismatch between database and BAM; "
                                    "expected size %i, found %i for contig %r"
                                    % (int(size), self._contigs[name]['Size'],
                                       name))

                row = {
                    'ID': name,
                    'Size': self._contigs[name]['Size'],
                    'Ns': self._contigs[name]['Ns'],
                    'Hits': hits,
                }

                handle.write('{ID}\t{Size}\t{Ns}\t{Hits}\n'.format(**row))

        CommandNode._setup(self, config, temp)


def hash_params(*args, **kwargs):
    return hashlib.md5(repr([args, kwargs])).hexdigest()
