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
import copy
import os

import pysam

import paleomix
import paleomix.resources

from paleomix.node import Node

import paleomix.common.fileutils as fileutils

import paleomix.tools.zonkey.parts.admixture as admixture
import paleomix.tools.zonkey.parts.nuclear as nuclear

from paleomix.tools.zonkey.common import \
    RSCRIPT_VERSION, \
    read_summary

from paleomix.nodes.samtools import SAMTOOLS_VERSION
from paleomix.nodes.raxml import RAXML_VERSION


class ReportNode(Node):
    def __init__(self, config, root, has_nuc, has_mt, dependencies=()):
        """

        Arguments:
          config -- Config object generated using paleomix.tools.zonkey.config.
          root -- Root folder containing current analysis.
          has_nuc -- True if a nuclear BAM was provided.
          has_mt -- True if a mitochondrial BAM was provided.
          dependencies -- Nodes for ReportNode to depend on.
        """

        self._root = root
        self._data = copy.deepcopy(config.database)
        self._report = AnalysisReport(config, root, has_nuc, has_mt)
        self._has_nuc = bool(has_nuc)
        self._has_mt = bool(has_mt)

        self._treemix_outgroup = config.treemix_outgroup
        self._treemix_k = config.treemix_k
        if self._treemix_k is None:
            self._treemix_k = '<automatic>'

        Node.__init__(self,
                      description="<Report -> %r>"
                      % (os.path.join(self._root, "report.html"),),
                      input_files=self._report.input_files(),
                      output_files=(os.path.join(self._root, "report.html"),
                                    os.path.join(self._root, "report.css")),
                      dependencies=dependencies)

    def _run(self, _config, temp):
        with open(os.path.join(temp, "report.html"), "w") as output_handle:
            revision = self._data.settings['Revision']
            header = _HTML_HEADER.format(Version=paleomix.__version__,
                                         Database=revision,
                                         Sidebar=self._build_sidebar())
            output_handle.write(header)

            self._write_intro_and_overview(output_handle)
            self._write_sample_description(output_handle)

            if self._has_nuc:
                self._write_admixture_estimates(output_handle)
                self._write_pca_plots(output_handle)
                self._write_treemix_plots(output_handle)

            if self._has_mt:
                self._write_mitochondrial_plots(output_handle)

            self._write_references(output_handle)

            output_handle.write(_HTML_FOOTER)

    def _teardown(self, config, temp):
        fileutils.make_dirs(self._root)

        fileutils.move_file(os.path.join(temp, "report.html"),
                            os.path.join(self._root, "report.html"))

        css_path = paleomix.resources.report("zonkey", "report.css")
        fileutils.copy_file(css_path, os.path.join(self._root, "report.css"))

    def _write_intro_and_overview(self, output_handle):
        output_handle.write(_SECTION_HEADER.format(name="intro",
                                                   title="Introduction"))
        output_handle.write(_INTRODUCTION)

        output_handle.write(_SECTION_HEADER.format(name="overview",
                                                   title="Analysis overview"))

        revision = self._data.settings['Revision']
        overview = _OVERVIEW_HEADER.format(DATABASE=revision,
                                           PYSAM=pysam.__version__,
                                           SAMTOOLS=_fmt_v(SAMTOOLS_VERSION),
                                           PLINK=_fmt_v(nuclear.PLINK_VERSION),
                                           RSCRIPT=_fmt_v(RSCRIPT_VERSION))
        output_handle.write(overview)

        if self._has_nuc:
            summary = self._report.snp_summary()
            output_handle.write(_OVERVIEW_NUCLEAR % summary)

        if self._has_mt:
            summary = self._report.mito_summary()
            output_handle.write(_OVERVIEW_MITOCHONDRIA % summary)

        if self._has_nuc:
            output_handle.write(_OVERVIEW_NUCLEAR_COVERAGE)

        output_handle.write(_OVERVIEW_FOOTER)

    def _write_sample_description(self, output_handle):
        output_handle.write(_SECTION_HEADER.format(name="samples",
                                                   title="Reference Panel"))

        output_handle.write(_SAMPLE_LIST_HEADER)

        last_group_2 = None
        last_group_3 = None
        for row in sorted(self._data.samples.itervalues(),
                          key=lambda row: (row["Group(2)"],
                                           row["Group(3)"],
                                           row["ID"])):

            row = dict(row)
            if last_group_2 != row["Group(2)"]:
                last_group_2 = row["Group(2)"]
                last_group_3 = row["Group(3)"]
            else:
                row["Group(2)"] = ""

                if last_group_3 != row["Group(3)"]:
                    last_group_3 = row["Group(3)"]
                else:
                    row["Group(3)"] = ""

            pub = row["Publication"]
            if pub.startswith("http"):
                row["Publication"] \
                    = '<a href="{0}">Link</a>'.format(pub.strip())
            elif row["Publication"].startswith("doi:"):
                pub = pub[4:].strip()
                url = "https://doi.org/{}".format(pub)
                row["Publication"] \
                    = 'doi:<a href="{0}">{1}</a>'.format(url, pub)

            output_handle.write(_SAMPLE_LIST_ROW.format(**row))

        output_handle.write("      </table>\n")

    def _write_admixture_estimates(self, output_handle):
        header = _SECTION_HEADER.format(name="admixture",
                                        title="Admixture Estimates")
        output_handle.write(header)

        admixture_v = _fmt_v(nuclear.ADMIXTURE_VERSION)
        overview = _ADMIXTURE_OVERVIEW.format(ADMIXTURE=admixture_v)
        output_handle.write(overview)

        for k_groups in (2, 3):
            summary_incl = self._build_admixture_cell(k_groups, True)
            summary_excl = self._build_admixture_cell(k_groups, False)

            output_handle.write(_ADMIXTURE_ROW.format(K=k_groups,
                                                      Incl_TS=summary_incl,
                                                      Excl_TS=summary_excl))

    def _build_admixture_cell(self, k_groups, incl_ts,
                              cutoff=admixture.CUTOFF):
        try:
            groups = self._report.admixture_results(k_groups, incl_ts)
        except admixture.AdmixtureError, error:
            return "<strong>ERROR: {}</strong".format(error)

        if sum((value >= cutoff) for _, value in groups) < 2:
            return "<strong>No admixture detected.</strong>"

        lines = [
            "<strong>Possible hybridization detected:</strong>",
            "<ul>",
        ]

        for group, value in groups:
            if value >= cutoff:
                name = " / ".join(sorted(group))

                lines.append("  <li>%s (%.2f%%)</li>" % (name, value * 100))

        lines.append("</ul>")

        if sum((value >= cutoff) for _, value in groups) != 2:
            return "\n            ".join(lines)

        percentiles = self._report.admixture_percentiles(data=self._data,
                                                         k_groups=k_groups,
                                                         incl_ts_k=incl_ts)

        if not ('Lower' in percentiles or 'Upper' in percentiles):
            lines.append("WARNING: Could not determine percentiles.")

            return "\n            ".join(lines)

        if 'Lower' not in percentiles:
            percentiles['Lower'] = percentiles['Upper']
            finale = \
                '; note that this is more simulated reads than what was ' \
                'processed in this analyses, potentially resulting in an ' \
                'overrestimating of percentages.'
        elif 'Upper' not in percentiles:
            percentiles['Upper'] = percentiles['Lower']
            finale = \
                '; note that this is fewer simulated reads than what was ' \
                'processed in this analyses, potentially resulting in an ' \
                'underrestimating of percentages'
        else:
            finale = '.'

        lower_pct = "%.1f" % ((1.0 - max(percentiles['Lower']['Upper'],
                                         percentiles['Upper']['Lower'])) *
                              100.0,)
        upper_pct = "%.1f" % ((1.0 - min(percentiles['Lower']['Upper'],
                                         percentiles['Upper']['Lower'])) *
                              100.0,)

        pct_range = lower_pct
        if lower_pct != upper_pct:
            pct_range = '%s - %s' % (lower_pct, upper_pct)

        lower_reads = min(percentiles['Lower']['NReads'],
                          percentiles['Upper']['NReads'])
        upper_reads = max(percentiles['Lower']['NReads'],
                          percentiles['Upper']['NReads'])

        reads_range = lower_reads
        if lower_reads != upper_reads:
            reads_range = '%s to %s' % (lower_reads, upper_reads)

        lines.append('Admixture results fall within %s percent of those '
                     'observed for simulated F1 %s / %s hybrids, based on '
                     '%s randomly selected reads%s'
                     % (pct_range,
                        percentiles['Sample1'],
                        percentiles['Sample2'],
                        reads_range, finale))

        return "\n            ".join(lines)

    def _write_pca_plots(self, output_handle):
        output_handle.write(_SECTION_HEADER.format(name="pca",
                                                   title="PCA Plots"))

        smartpca_v = _fmt_v(nuclear.SMARTPCA_VERSION)
        output_handle.write(_PCA_SECTION.format(SMARTPCA=smartpca_v))

    def _write_treemix_plots(self, output_handle):
        output_handle.write(_SECTION_HEADER.format(name="treemix",
                                                   title="Treemix Plots"))

        outgroups = ""
        if self._treemix_outgroup:
            outgroups = ", ".join(map(repr, self._treemix_outgroup))
            outgroups = ". The tree was rooted on the clade containing the " \
                        "sample(s) %s" % (outgroups)

        treemix_v = _fmt_v(nuclear.TREEMIX_VERSION)
        overview = _TREEMIX_OVERVIEW.format(treemix_k=self._treemix_k,
                                            treemix_outgroup=outgroups,
                                            TREEMIX=treemix_v)
        output_handle.write(overview)

        for prefix in ("incl_ts", "excl_ts"):
            output_handle.write("<h2>%s</h2>\n" % (_TS_LABELS[prefix],))

            for n_edges in (0, 1):
                variance_file = os.path.join(self._root,
                                             "figures",
                                             "treemix",
                                             "%s_%i_variance.txt")

                with open(variance_file % (prefix, n_edges)) as handle:
                    variance = handle.read().strip()

                treemix_row = _TREEMIX_TREE_ROW.format(Prefix=prefix,
                                                       Edges=n_edges,
                                                       Variance=variance)
                output_handle.write(treemix_row)

    def _write_mitochondrial_plots(self, output_handle):
        header = _SECTION_HEADER.format(name="mitochondria",
                                        title="Mitochondrial Phylogeny")
        output_handle.write(header)
        raxml_v = _fmt_v(RAXML_VERSION)
        output_handle.write(_MITOCONDRIAL_SECTION.format(RAXML=raxml_v))

    def _write_references(self, output_handle):
        header = _SECTION_HEADER.format(name="references",
                                        title="References")
        output_handle.write(header)
        output_handle.write(_REFERENCES)

    def _build_sidebar(self):
        lines = [_SIDEBAR_HEADER]

        if self._has_nuc:
            lines.append(_SIDEBAR_NUCLEAR)

        if self._has_mt:
            lines.append(_SIDEBAR_MITO)

        lines.append(_SIDEBAR_FOOTER)

        return "\n".join(lines)


class AnalysisReport(object):
    def __init__(self, config, root, has_nuc, has_mt):
        self._has_nuc = bool(has_nuc)
        self._has_mt = bool(has_mt)
        self._filtered = bool(config.indep)
        self._config = config
        self._root = root
        self._data = config.database

    def input_files(self):
        input_files = [self._config.tablefile]
        if self._has_nuc:
            # Summary file generated while building TPED files
            input_files.append(os.path.join(self._root,
                                            "results",
                                            "plink",
                                            "common.summary"))

            for postfix in ('incl_ts', 'excl_ts'):
                admix_root = os.path.join(self._root, "results", "admixture")

                if self._filtered:
                    # Required to count number of SNPs included after filtering
                    input_files.append(os.path.join(self._root,
                                                    "results",
                                                    "plink",
                                                    postfix + ".bim"))

                # Include files showing proproation of ancestral populations,
                # which are required to build admixture figures in the reports.
                for k_groups in (2, 3):
                    input_files.append(os.path.join(admix_root,
                                                    "%s.%i.Q" % (postfix,
                                                                 k_groups)))

                # Include files tabulating variance explained by models
                figures_path = os.path.join(self._root, "figures", "treemix")
                for n_edges in (0, 1):
                    variance_path = os.path.join(figures_path,
                                                 "%s_%i_variance.txt"
                                                 % (postfix, n_edges))

                    input_files.append(variance_path)

        if self._has_mt:
            input_files.append(os.path.join(self._root,
                                            "results",
                                            "mitochondria",
                                            "sequences.summary"))

        return input_files

    def snp_summary(self):
        summary = read_summary(os.path.join(self._root,
                                            "results",
                                            "plink",
                                            "common.summary"))

        if self._filtered:
            for postfix in ('incl_ts', 'excl_ts'):
                key = "n_sites_%s" % (postfix,)

                filename = os.path.join(self._root, "plink", postfix + ".bim")
                with open(filename) as handle:
                    n_used = sum(1 for _ in handle)

                summary[key] = "%i of %i" % (n_used, summary[key])

        return summary

    def mito_summary(self):
        return read_summary(os.path.join(self._root,
                                         "results",
                                         "mitochondria",
                                         "sequences.summary"))

    def admixture_results(self, k_groups, incl_ts,
                          cutoff=admixture.CUTOFF):
        prefix = "incl_ts" if incl_ts else "excl_ts"
        filename = os.path.join(self._root, "results", "admixture",
                                "%s.%i.Q" % (prefix, k_groups))

        return admixture.read_admixture_results(filename=filename,
                                                data=self._config.database,
                                                k_groups=k_groups,
                                                cutoff=cutoff)

    def admixture_percentiles(self, data, k_groups, incl_ts_k,
                              cutoff=admixture.CUTOFF):
        try:
            results = self.admixture_results(k_groups, incl_ts_k, cutoff)
        except admixture.AdmixtureError:
            return

        groups = [group for group, value in results if value >= cutoff]
        if len(groups) != 2:
            return

        sample1, = groups[0]
        sample2, = groups[1]
        delta = abs(max(value for _, value in results) - 0.5)
        summary = read_summary(os.path.join(self._root,
                                            "results",
                                            "plink",
                                            "common.summary"))

        return admixture.get_percentiles(data=data,
                                         sample1=sample1,
                                         sample2=sample2,
                                         nreads=summary['n_reads_used'],
                                         k_groups=k_groups,
                                         has_ts=incl_ts_k,
                                         value=delta)


def _fmt_v(requirement):
    return ".".join(map(str, requirement.version))


###############################################################################

_TS_LABELS = {
    "incl_ts": "Including transitions",
    "excl_ts": "Excluding transitions",
}


###############################################################################

_HTML_HEADER = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>PALEOMIX Zonkey v{Version} - db rev. {Database}</title>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<link href="report.css" rel="stylesheet" type="text/css" />
</head>
<body>
<a name="top" id="top"></a>
<center>
  <div id="header">
    <h1>PALEOMIX Zonkey v{Version} - db rev. {Database}</h1>
    <h2>A pipeline for detection of F1 hybrids in equids.</h2>
  </div>
  <div id="content">
    <p class="introduction">
      Schubert M, Ermini L, Sarkissian CD, J&oacute;nsson H, Ginolhac A,
      Schaefer R, Martin MD, Fern&aacute;ndez R, Kircher M, McCue M,
      Willerslev E, and Orlando L. "<strong>Characterization of ancient and
      modern genomes by SNP detection and phylogenomic and metagenomic analysis
      using PALEOMIX</strong>". Nat Protoc. 2014 May;9(5):1056-82. doi:<a
      href="https://doi.org/10.1038/nprot.2014.063">10.1038/nprot.2014.063</a>.
      Epub 2014 Apr 10. PubMed PMID: <a
      href="http://www.ncbi.nlm.nih.gov/pubmed/24722405">24722405</a>.
    </p>
{Sidebar}
    <div id="mainbar">
"""

_SECTION_HEADER = """      <h1><a name="{name}" id="{name}"></a>{title}</h1>
"""

_INTRODUCTION = """
      <div>
        <div>
          The Zonkey Pipeline is a easy-to-use pipeline designed for the
          analyses of low-coverage, ancient DNA derived from historical
          equid samples, with the purpose of determining the species of
          the sample, as well as determining possible hybridization between
          horses, zebras, and asses. This is accomplished by comparing one
          or more samples aligned against the <em>Equus caballus</em> 2.0
          reference sequence with a reference panel of modern equids,
          including wild and domesticated equids.
        </div>
        <br/>
        <div>
          For more information, please refer to the
          <a href="http://paleomix.readthedocs.org/en/latest/zonkey_pipeline/index.html">
            the documentation for the Zonkey pipeline
          </a>
          or
          <a href="http://paleomix.readthedocs.org/en/latest/">
            the documentation for the PALEOMIX pipeline,
          </a>
          on which the Zonkey pipeline is based.
        </div>
        <br/>
"""

_OVERVIEW_HEADER = """
        <div>
          Zonkey run using database rev. {DATABASE}. Data processed using
          <a href="https://github.com/pysam-developers/pysam">pysam</a> v{PYSAM},
          <a href="https://samtools.github.io/">SAMTools</a> v{SAMTOOLS}
          [<em>Li <em>et al.</em> 2009</em>] and
          <a href="http://pngu.mgh.harvard.edu/purcell/plink/">PLINK</a> v{PLINK}
          [<em>Purcell <em>et al.</em> 2007</em>]; plotting was carried out using
          <a href="https://www.r-project.org/">R</a> v{RSCRIPT}. Additional
          tools listed below.
        </div>
        <br/>
        <div style="display:table;width:100%;">
          <div style="display:table-cell;">
"""

_OVERVIEW_NUCLEAR = """
            <strong>Nuclear report from '<em>%(filename)s</em>'</strong>

            <table style="width:95%%">
              <tr>
                <td style="width:50%%;">Number of reads processed:</td>
                <td>%(n_reads)s</td>
              </tr>
              <tr>
                <td>Number of reads overlapping SNPs:</td>
                <td>%(n_reads_used)s</td>
              </tr>
              <tr>
                <td>Number of SNPs used (incl. transitions):</td>
                <td>%(n_sites_incl_ts)s</td>
              </tr>
              <tr>
                <td>Number of SNPs used (excl. transitions):</td>
                <td>%(n_sites_excl_ts)s</td>
              </tr>
            </table>
"""

_OVERVIEW_MITOCHONDRIA = """
            <br>

            <h4>Mitochondrial report from '<em>%(filename)s</em>'</h4>

            <table style="width:95%%">
              <tr>
                <td style="width:50%%;">Reference sequence used:</td>
                <td>%(sequence_name)s</td>
              </tr>
              <tr>
                <td>Reference sequence length:</td>
                <td>%(sequence_len)s</td>
              </tr>
              <tr>
                <td>Number of sites covered:</td>
                <td>%(covered_sites)s</td>
              </tr>
              <tr>
                <td>Percentage of sites covered:</td>
                <td>%(covered_pct)s</td>
              </tr>
              <tr>
                <td>Mean coverage per site:</td>
                <td>%(mean_coverage)s</td>
              </tr>
            </table>
"""


_OVERVIEW_NUCLEAR_COVERAGE = """
          </div>
          <div style="display:table-cell;width:25%;max-width:350px;">
           <div>
             <strong>Autosomes vs. sex-chromosomes:</strong>
           </div>
           <div>
             <a href="figures/coverage/coverage.pdf">
               <img src="figures/coverage/coverage.png"
                    style="vertical-align:top;horizontal-align:center;">
             </a>
           </div>
"""

_OVERVIEW_FOOTER = """
          </div>
      </div>
"""

_SAMPLE_LIST_HEADER = """
      <table summary="List of samples in the reference panel.">
        <tr>
          <th>Group(2)</th>
          <th>Group(3)</th>
          <th>ID</th>
          <th>Species</th>
          <th>Sex</th>
          <th>Sample Name</th>
          <th>Publication</th>
        </tr>
"""

_SAMPLE_LIST_ROW = """
        <tr>
          <td>{Group(2)}
          <td>{Group(3)}
          <td>{ID}</td>
          <td><em>{Species}</em></td>
          <td>{Sex}</th>
          <td>{SampleID}</td>
          <td>{Publication}</td>
        </tr>
"""


_ADMIXTURE_OVERVIEW = """
      <p>
        Admixture proportions estimated using
        <a href="https://www.genetics.ucla.edu/software/admixture/">ADMIXTURE</a>
        v{ADMIXTURE} <em>[Alexander <em>et al.</em> 2009]</em>, using default
        parameters.
      </p>
"""


_ADMIXTURE_ROW = """
      <table summary="Admixture between sample and reference populations.">
        <tr>
          <th>{K} ancestral groups</th>
          <th>{K} ancestral groups, excluding transitions</th>
        </tr>
        <tr>
          <td>
            <a href="figures/admixture/incl_ts_k{K}.pdf">
                <img style="width:95%" src="figures/admixture/incl_ts_k{K}.png"
                     alt="Admixture for k={K}, including transitions." />
            </a>
            <br/>
            <br/>
            {Incl_TS}
          </td>
          <td>
            <a href="figures/admixture/excl_ts_k{K}.pdf">
                <img style="width:95%" src="figures/admixture/excl_ts_k{K}.png"
                     alt="Admixture for k={K}, excluding transitions." />
            </a>
            <br/>
            <br/>
            {Excl_TS}
          </td>
        </tr>
      </table>
"""


_PCA_SECTION = """
      <p>
        Principal Component Analysis carried out using SmartPCA v{SMARTPCA},
        from the <a
        href="http://www.hsph.harvard.edu/alkes-price/software/">EIGENSOFT</a>
        toolkit.
      </p>

      <table summary="PCA plots comparing sample with the reference panel.">
        <tr>
          <th>Including transitions</th>
          <th>Excluding transitions</th>
        </tr>
        <tr>
          <td>
            <a href="figures/pca/incl_ts.pdf">
                <img style="width:95%" src="figures/pca/incl_ts.png"
                     alt="PCA plot, including transitions." />
            </a>
          </td>
          <td>
            <a href="figures/pca/excl_ts.pdf">
                <img style="width:95%" src="figures/pca/excl_ts.png"
                     alt="PCA plot, excluding transitions." />
            </a>
          </td>
        </tr>
      </table>
"""


_TREEMIX_OVERVIEW = """
      <p>
        Detection of population mixture using
        <a href="https://bitbucket.org/nygcresearch/treemix/wiki/Home">TreeMix</a>
        v{TREEMIX} <em>[Pickrell and Pritchard 2012]</em>; parameters were -k
        {treemix_k}; -global; and supervised estimation using ancestral groups
        listed in the Reference Panel{treemix_outgroup}.
      </p>
"""


_TREEMIX_TREE_ROW = """
      <table summary="Treemix plots, for {Edges} edge(s).">
        <tr>
          <th>Edges = {Edges}</th>
          <th>Residuals</th>
        </tr>
        <tr>
          <td>
            <a href="figures/treemix/{Prefix}_{Edges}_tree.pdf">
                <img style="width:95%"
                     src="figures/treemix/{Prefix}_{Edges}_tree.png"
                     alt="Treemix plot, {Edges} edge(s), incl. transitions." />
            </a>
          </td>
          <td>
            <a href="figures/treemix/{Prefix}_{Edges}_residuals.pdf">
                <img style="width:95%"
                     src="figures/treemix/{Prefix}_{Edges}_residuals.png"
                     alt="Treemix plot, {Edges} edge(s), excl. transitions." />
            </a>
          </td>
        </tr>
      </table>

      <p>
        Variance explained by model = {Variance}.
      </p>
"""


_MITOCONDRIAL_SECTION = """
      <p>
        Phylogenetic inference performed using RAxML
        v{RAXML} [<em>Stamatakis 2006</em>].
      <p>

      <div style="text-align:center;">
        <a href="figures/mitochondria/mito_phylo.pdf">
          <img src="figures/mitochondria/mito_phylo.png"
               alt="Mitochondrial maximum likelihood phylogeny." />
        </a>
      </div>
"""


_REFERENCES = """
        <p>
          <ul>
            <li>
              Alexander <em>et al</em>. "<strong>Fast model-based estimation of
              ancestry in unrelated individuals</strong>". <em>Genome Res</em>.
              2009 Sep;19(9):1655-64.
                doi:<a href="https://doi.org/10.1101/gr.094052.109">
                    10.1101/gr.094052.109</a>.
                PMID:<a href="http://www.ncbi.nlm.nih.gov/pubmed/19648217">
                  19648217</a>.
            </li>
            <li>
              Li <em>et al</em>. "<strong>The Sequence Alignment/Map format and
              SAMtools</strong>. <em>Bioinformatics</em>. 2009 Aug
              15;25(16):2078-9.
                doi:<a href="https://doi.org/10.1093/bioinformatics/btp352">
                    10.1093/bioinformatics/btp352</a>.
                PMID:<a href="http://www.ncbi.nlm.nih.gov/pubmed/19505943">
                  19505943</a>.
            </li>
            <li>
              Pickrell and Pritchard. "<strong>Inference of population splits
              and mixtures from genome-wide allele frequency data</strong>".
              <em>PLoS Genet</em>. 2012;8(11):e1002967.
              doi:<a href="https://doi.org/10.1371/journal.pgen.1002967">
                  10.1371/journal.pgen.1002967</a>.
              PMID:<a href="http://www.ncbi.nlm.nih.gov/pubmed/23166502">
                23166502</a>.
            </li>
            <li>
              Purcell <em>et al</em>. "<strong>PLINK: a tool set for whole-
              genome association and population-based linkage
              analyses</strong>". <em>Am J Hum Genet</em>. 2007
              Sep;81(3):559-75.
              PMID:<a href="http://www.ncbi.nlm.nih.gov/pubmed/17701901">
                17701901</a>.
            </li>
            <li>
              Stamatakis. "<strong>RAxML-VI-HPC: maximum likelihood-based
              phylogenetic analyses with thousands of taxa and mixed
              models</strong>". <em>Bioinformatics</em>. 2006 Nov
              1;22(21):2688-90. Epub 2006 Aug 23.
              doi:<a href="https://doi.org/10.1093/bioinformatics/btl446">
                  10.1093/bioinformatics/btl446</a>.
              PMID:<a href="http://www.ncbi.nlm.nih.gov/pubmed/16928733">
                16928733</a>.
          </ul>
        </p>
"""


_HTML_FOOTER = """
    </div>
  </div>
  <div id="footer">
    This report is based on the PLAIN 1.0 design by
    <a href="http://www.sixshootermedia.com/">6ix Shooter Media</a>,
    Creative Commons license.<br/>
  </div>
</center>
</body>
</html>
"""

###############################################################################

_SIDEBAR_HEADER = """
    <div id="sidebar">
      <h1>Contents</h1>
      <div class="submenu">
        <a href="#">Top</a>
        <a href="#intro">Introduction</a>
        <a href="#overview">Analysis overview</a>
        <a href="#samples">Reference Panel</a>
"""

_SIDEBAR_NUCLEAR = """
        <a href="#admixture">Admixture Estimates</a>
        <a href="#pca">PCA Plots</a>
        <a href="#treemix">Treemix Analyses</a>
"""

_SIDEBAR_MITO = """
        <a href="#mito_phylo">MT Phylogeny</a>
"""

_SIDEBAR_FOOTER = """
        <a href="#references">References</a>
      </div>
    </div>
"""
