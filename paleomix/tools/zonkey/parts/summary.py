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
import os

import paleomix
import paleomix.resources

from paleomix.node import Node

import paleomix.common.fileutils as fileutils

import paleomix.tools.zonkey.parts.admixture as admixture

from paleomix.tools.zonkey.parts.report import AnalysisReport


class SummaryNode(Node):
    def __init__(self, config, dependencies=()):
        self._root = config.destination
        self._data = config.database
        self._samples = config.samples
        self._sample_keys = self._samples.keys()

        input_files = set()
        self._reports = {}
        for sample, info in self._samples.iteritems():
            report = AnalysisReport(config=config,
                                    root=os.path.join(self._root, sample),
                                    has_nuc="Nuc" in info["Files"],
                                    has_mt="Mito" in info["Files"])

            input_files.update(report.input_files())
            self._reports[sample] = report

        output_prefix = os.path.join(self._root, "summary")
        Node.__init__(self,
                      description="<SummaryReport -> %r>"
                      % (output_prefix + '.html',),
                      input_files=input_files,
                      output_files=(output_prefix + '.html',
                                    output_prefix + '.css'),
                      dependencies=dependencies)

    def _run(self, _config, temp):
        with open(os.path.join(temp, "summary.html"), "w") as output_handle:
            menu_entries = self._build_sidemenu()
            html_header = _HTML_HEADER.format(Version=paleomix.__version__,
                                              MenuEntries=menu_entries)
            output_handle.write(html_header)

            for sample in sorted(self._samples):
                self._write_sample_overview(output_handle, sample)

            output_handle.write(_HTML_FOOTER)

    def _teardown(self, config, temp):
        fileutils.make_dirs(self._root)

        fileutils.move_file(os.path.join(temp, "summary.html"),
                            os.path.join(self._root, "summary.html"))

        css_path = paleomix.resources.report("zonkey", "report.css")
        fileutils.copy_file(css_path, os.path.join(self._root, "summary.css"))

    def _build_sidemenu(self):
        lines = []
        for sample in sorted(self._samples):
            lines.append('        <a href="#sample_{0}">{0}</a>'.format(sample))
        return "\n".join(lines)

    def _write_sample_overview(self, handle, sample):
        info = self._samples[sample]
        sample_idx = self._sample_keys.index(sample)

        admixture = "&nbsp;"
        if "Nuc" in info["Files"]:
            admixture = self._read_admixture_results(sample)

        handle.write(_SAMPLE_HEADER.format(sample=sample,
                                           admixture=admixture))

        handle.write("""    <ul id="tabs">""")
        handle.write(_SAMPLE_TAB_SELECTED.format(sample_idx=sample_idx,
                                                 page=1, title="Overview"))
        handle.write(_SAMPLE_TAB.format(sample_idx=sample_idx,
                                        page=2, title="Sample Data"))

        if "Nuc" in info["Files"]:
            handle.write(_SAMPLE_TAB.format(sample_idx=sample_idx,
                                            page=3, title="PCA Plots"))
            handle.write(_SAMPLE_TAB.format(sample_idx=sample_idx,
                                            page=4, title="Treemix Plots"))

        if "Mito" in info["Files"]:
            handle.write(_SAMPLE_TAB.format(sample_idx=sample_idx,
                                            page=5,
                                            title="Mitochondrial Phylogeny"))

        handle.write("""    </ul>""")

        self._write_overview(handle, sample, sample_idx)

        if "Nuc" in info["Files"]:
            handle.write(_SAMPLE_PAGE_PCA.format(sample=sample,
                                                 sample_idx=sample_idx))
            handle.write(_SAMPLE_PAGE_TREEMIX.format(sample=sample,
                                                     sample_idx=sample_idx))

        if "Mito" in info["Files"]:
            handle.write(_SAMPLE_PAGE_MITO_PHYLO.format(sample=sample,
                                                        sample_idx=sample_idx))

        handle.write(_SAMPLE_FOOTER)

    def _write_overview(self, output_handle, sample, sample_idx,
                        cutoff=admixture.CUTOFF):
        info = self._samples[sample]
        report = self._reports[sample]

        n_tests = 4
        n_pos_tests = 0

        if "Nuc" in info["Files"]:

            for postfix in ('incl_ts', 'excl_ts'):
                admix_root = os.path.join(self._root, sample,
                                          "results", "admixture")

                for k_groups in (2, 3):
                    filename = os.path.join(admix_root, "%s.%i.Q" % (postfix,
                                                                     k_groups))

                    result = admixture.read_admixture_results(filename,
                                                              self._data,
                                                              k_groups)

                    if sum(value >= cutoff for _, value in result) > 1:
                        n_pos_tests += 1

            if n_pos_tests:
                tmpl = _SAMPLE_OVERVIEW_INCL_NUCL_POSITIVE
            else:
                tmpl = _SAMPLE_OVERVIEW_INCL_NUCL_NEGATIVE
        else:
            tmpl = _SAMPLE_OVERVIEW_EXCL_NUCL

        output_handle.write(tmpl.format(sample_idx=sample_idx,
                                        n_tests=n_tests,
                                        n_pos_tests=n_pos_tests))
        output_handle.write(_SAMPLE_OVERVIEW_FOOTER.format(sample=sample))

        output_handle.write(_SAMPLE_DATA_HEADER.format(sample_idx=sample_idx))

        if "Nuc" in info["Files"]:
            summary = report.snp_summary()
            output_handle.write(_SAMPLE_DATA_NUCLEAR % summary)

        if "Mito" in info["Files"]:
            summary = report.mito_summary()
            output_handle.write(_SAMPLE_DATA_MITOCHONDRIA % summary)

        output_handle.write(_SAMPLE_DATA_FOOTER)

    def _read_admixture_results(self, sample, cutoff=admixture.CUTOFF):
        lines = []
        lines.append('            <table summary="Admixture overview for sample {}" style="width:125px;">'.format(sample.replace('"', "")))
        lines.append("              <tr>")
        for postfix in ('incl_ts', 'excl_ts'):
            admix_root = os.path.join(self._root, sample,
                                      "results", "admixture")

            for k_groups in (2, 3):
                filename = os.path.join(admix_root, "%s.%i.Q" % (postfix, k_groups))

                lines.append("                <td>")
                try:
                    ancestral_groups = admixture.read_admixture_results(filename, self._data, k_groups, cutoff)
                    lines.extend(self._build_admixture_figure([value for _, value in ancestral_groups]))
                except admixture.AdmixtureError:
                    lines.append('                  <div style="height:100px;background:gray"></div>')
                lines.append("                </td>")
        lines.append("              </tr>")
        lines.append("            </table>\n")

        return "\n".join(lines)

    @classmethod
    def _build_admixture_figure(cls, fractions, max_height=100):
        lines = []
        if len(fractions) == 2:
            colors = ("red", "green")
        elif len(fractions) == 3:
            colors = ("red", "green", "blue")
        else:
            raise RuntimeError("Unexpected number of fractions: %r"
                               % (fractions,))

        values = [value for value in fractions if value >= admixture.CUTOFF]

        tmpl = ' ' * 18 + '<div style="height:{}px;background:{}"></div>'
        for index, value in enumerate(values):
            height = round((float(value) / sum(values)) * max_height)
            lines.append(tmpl.format(height, colors[index]))

        return lines


###############################################################################

_TS_LABELS = {
    "incl_ts": "Including transitions",
    "excl_ts": "Excluding transitions",
}


###############################################################################

_HTML_HEADER = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
    <title>PALEOMIX Zonkey v{Version}</title>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
    <link href="summary.css" rel="stylesheet" type="text/css" />

    <script type="text/javascript">
        function selectTab(sample, page) {{
            for (var i = 1; i <= 5; ++i) {{
                elem_id = "sample" + sample.toString() + ".link" + i.toString()

                var elem = document.getElementById(elem_id);
                if (elem) {{
                    elem.className = (i == page) ? 'selected' : '';
                }}

                elem_id = "sample" + sample.toString() + ".page" + i.toString()
                var elem = document.getElementById(elem_id);
                if (elem) {{
                    if (i == page) {{
                        if (page == 1) {{
                            elem.className = 'tabContent small';
                        }} else {{
                            elem.className = 'tabContent';
                        }}
                    }} else {{
                        elem.className = 'tabContent hide';
                    }}
                }}
            }}
        }}
    </script>
</head>
<body>
  <a name="top" id="top"></a>
  <center>
    <div id="header">
      <h1>PALEOMIX Zonkey {Version}</h1>
      <h2>A pipeline for detection of F1 hybrids in equids.</h2>
    </div>

    <div id="content" style="width:1050px;">
      <p class="introduction">
        Schubert M, Ermini L, Sarkissian CD, J&oacute;nsson H, Ginolhac A,
        Schaefer R, Martin MD, Fern&aacute;ndez R, Kircher M, McCue M,
        Willerslev E, and Orlando L. "<strong>Characterization of ancient and
        modern genomes by SNP detection and phylogenomic and metagenomic
        analysis using PALEOMIX</strong>". Nat Protoc. 2014 May;9(5):1056-82.
        doi:<a href="https://doi.org/10.1038/nprot.2014.063">
            10.1038/nprot.2014.063
        </a>.
        Epub 2014 Apr 10. PubMed PMID:
        <a href="http://www.ncbi.nlm.nih.gov/pubmed/24722405">24722405</a>.
      </p>
      <div id="sidebar">
        <h1>Contents</h1>
        <div class="submenu">
          <a href="#">Top</a>
          <a href="#intro">Introduction</a>
          {MenuEntries}
        </div>
      </div>
      <div id="mainbar">
        <h1><a name="introduction" id="introduction"></a>Introduction</h1>
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
        </div>
        <h1><a name="samples" id="samples"></a>Sample Overview</h1>
"""


_SAMPLE_HEADER = """
    <h2 style="margin: 0px; padding: 0px;">
      <a name="sample_{sample}" id="sample_{sample}"
         href="{sample}/report.html">{sample}</a>
    </h2>

    <table style="margin: 0px; padding: 0px;">
      <tr>
        <td style="width:150px; vertical-align:top; padding-top:20px;">
{admixture}
        </td>
        <td>
"""

_SAMPLE_FOOTER = """
        </td>
      </tr>
    </table>
"""

_SAMPLE_TAB_SELECTED = """
      <li>
        <a id="sample{sample_idx}.link{page}"
           href="javascript:selectTab({sample_idx}, {page});"
           class="selected">{title}</a>
      </li>
"""

_SAMPLE_TAB = """
      <li>
        <a id="sample{sample_idx}.link{page}"
           href="javascript:selectTab({sample_idx}, {page});">{title}</a>
      </li>
"""


_SAMPLE_OVERVIEW_INCL_NUCL_POSITIVE = """
    <div class="tabContent small" id="sample{sample_idx}.page1">
      <div>
        Possible admixture detected in {n_pos_tests} of {n_tests} tests.
      </div>
"""


_SAMPLE_OVERVIEW_INCL_NUCL_NEGATIVE = """
    <div class="tabContent small" id="sample{sample_idx}.page1">
      <div>
        No admixture detected in nuclear genome.
      </div>
"""


_SAMPLE_OVERVIEW_EXCL_NUCL = """
    <div class="tabContent small" id="sample{sample_idx}.page1">
      <div>
        Nuclear BAM not provided; no admixture were tests perfomed.
      </div>
"""


_SAMPLE_OVERVIEW_FOOTER = """
      <div style="text-align:center;">
          <a href="{sample}/report.html">(Click here for the full report)</a>
      </div>
    </div>
"""


_SAMPLE_DATA_HEADER = """
    <div class="tabContent hide" id="sample{sample_idx}.page2">
      <ul>
"""


_SAMPLE_DATA_NUCLEAR = """
      <strong>Nuclear report:</strong>
      <ul>
        <li>BAM file: <em>%(filename)s</em></li>
        <li>Number of SNPs used (incl. transitions):
            <em>%(n_sites_incl_ts)s</em></li>
        <li>Number of SNPs used (excl. transitions):
            <em>%(n_sites_excl_ts)s</em></li>
      </ul>

      <br>
"""


_SAMPLE_DATA_MITOCHONDRIA = """
      <strong>Mitochondrial report:</strong>
      <ul>
        <li>BAM file: <em>%(filename)s</em></li>
        <li>Percentage of sites covered: <em>%(covered_pct)s</em></li>
        <li>Mean coverage per site: <em>%(mean_coverage)s</em></li>
      </ul>
"""


_SAMPLE_DATA_FOOTER = """
      </ul>
    </div>
"""


_SAMPLE_PAGE_PCA = """
    <div class="tabContent hide" id="sample{sample_idx}.page3">
      <table>
        <tr style="background-color:#f1f0ee">
          <td style="text-align: center;">
            <div class="thumbnail">
              <strong>Including transitions</strong>

              <div class="image">
                <a href="{sample}/figures/pca/incl_ts.pdf">
                  <img style="width:13em"
                       src="{sample}/figures/pca/incl_ts.png"
                       alt="PCA plot for {sample}, including transitions."/>
                </a>
              </div>
            </div>
          </td>
          <td style="text-align: center;">
            <div class="thumbnail">
              <strong>Excluding transitions</strong>
              <div class="image">
                <a href="{sample}/figures/pca/excl_ts.pdf">
                  <img style="width:13em"
                       src="{sample}/figures/pca/excl_ts.png"
                       alt="PCA plot for {sample}, excluding transitions."/>
                </a>
              </div>
            </div>
          </td>
        </tr>
      </table>
    </div>
"""


_SAMPLE_PAGE_TREEMIX = """
    <div class="tabContent hide" id="sample{sample_idx}.page4">
      <table>
        <tr style="background-color:#f1f0ee">
          <td style="text-align: center;">
            <div class="thumbnail">
              <strong>Including transitions</strong>
              <div class="image">
                <a href="{sample}/figures/treemix/incl_ts_1_tree.pdf">
                  <img style="width:13em"
                       src="{sample}/figures/treemix/incl_ts_1_tree.png"
                       alt="Treemix plot for {sample}, including transitions,
                            with one migration edge."/>
                </a>
              </div>
            </div>
          </td>
          <td style="text-align: center;">
            <div class="thumbnail">
              <strong>Excluding transitions</strong>
              <div class="image">
                <a href="{sample}/figures/treemix/excl_ts_1_tree.pdf">
                  <img style="width:13em"
                       src="{sample}/figures/treemix/excl_ts_1_tree.png"
                       alt="Treemix plot for {sample}, excluding transitions,
                            with one migration edge."></img>
                </a>
              </div>
            </div>
          </td>
        </tr>
      </table>
    </div>
"""


_SAMPLE_PAGE_MITO_PHYLO = """
    <div class="tabContent hide" id="sample{sample_idx}.page5">
      <div class="thumbnail" style="text-align: center;">
        <div class="image">
          <a href="{sample}/figures/mitochondria/mito_phylo.pdf">
            <img style="width:15em"
                 src="{sample}/figures/mitochondria/mito_phylo.png"
                 alt="Mitochondrial phylogeny for {sample}."/>
          </a>
        </div>
      </div>
    </div>
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
