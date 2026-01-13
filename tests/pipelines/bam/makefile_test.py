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

import copy
import textwrap
from pathlib import Path
from typing import Any

import pytest

from paleomix.common.makefile import MakefileError
from paleomix.common.resources import read_template
from paleomix.pipelines.bam.makefile import read_makefiles


def _read_makefile(filepath: Path) -> object:
    (result,) = read_makefiles([filepath])

    return result


########################################################################################
# Test the validity of the makefile templates

# Options specified in a basic YAML file
TEMPLATE_BAM_OPTIONS = {
    "AdapterRemoval": {
        "Version": 3,
        "--collapse": None,
        "--minlength": 25,
    },
    "Aligners": {
        "Program": "BWA",
        "BWA": {
            "Algorithm": "mem",
            "FilterUnmappedReads": True,
            "MinQuality": 0,
            "UseSeed": True,
        },
        "Bowtie2": {
            "--very-sensitive": None,
            "FilterUnmappedReads": True,
            "MinQuality": 0,
        },
    },
    "ExcludeReads": {
        "Collapsed": False,
        "CollapsedTruncated": False,
        "Paired": False,
        "Single": False,
        "Singleton": False,
    },
    "Features": {
        "Coverage": True,
        "Depths": True,
        "PCRDuplicates": "filter",
        "Summary": True,
        "mapDamage": False,
    },
    "Platform": "ILLUMINA",
    "QualityOffset": 33,
    "mapDamage": {"--downsample": 100000},
}


TEMPLATE_BAM = {
    "Genomes": {
        "NAME_OF_PREFIX": {
            "Path": "PATH_TO_PREFIX.fasta",
            "Name": "NAME_OF_PREFIX",
            "IndexFormat": ".bai",
        },
    },
    "Options": copy.deepcopy(TEMPLATE_BAM_OPTIONS),
    "Samples": {
        "NAME_OF_SAMPLE": {
            "NAME_OF_LIBRARY": {
                "NAME_OF_LANE": [
                    {
                        "Path": ("PATH_WITH_WILDCARDS", None),
                        "Shortname": "PATH_WITH_WILDCARDS",
                        "Options": copy.deepcopy(TEMPLATE_BAM_OPTIONS),
                        "Type": "Untrimmed",
                    }
                ]
            }
        }
    },
}


def _trim_template() -> dict[str, Any]:
    data = copy.deepcopy(TEMPLATE_BAM)
    del data["Genomes"]["NAME_OF_PREFIX"]
    del data["Options"]["Aligners"]["Bowtie2"]["--very-sensitive"]
    del data["Options"]["mapDamage"]["--downsample"]

    for dd in data["Samples"]["NAME_OF_SAMPLE"]["NAME_OF_LIBRARY"]["NAME_OF_LANE"]:
        del dd["Options"]["Aligners"]["Bowtie2"]["--very-sensitive"]
        del dd["Options"]["mapDamage"]["--downsample"]

    return data


TEMPLATE_TRIM = _trim_template()


def test_basic__read_makefiles__trim_template(tmp_path: Path) -> None:
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        # Minimal template containing options required for trim pipeline
        handle.write(read_template("bam_head.yaml"))
        # Example samples (not used when auto generating a YAML file)
        handle.write(read_template("bam_samples.yaml"))

    expected = copy.deepcopy(TEMPLATE_TRIM)
    expected["Filename"] = filepath

    # The template options should always be valid
    assert _read_makefile(filepath) == expected


# The BAM template is split into parts to support both trimming and BAM pipelines
TEMPLATE_FILES = [
    # Minimal template containing options required for trim pipeline
    "bam_head.yaml",
    # Options used for the full BAM pipeline
    "bam_options.yaml",
    # Genome template (not used in trimming pipeline)
    "bam_prefixes.yaml",
    # Example samples (not used when auto generating a YAML file)
    "bam_samples.yaml",
]


def test_basic__read_makefiles__bam_template(tmp_path: Path) -> None:
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES:
            handle.write(read_template(filename))

    expected = copy.deepcopy(TEMPLATE_BAM)
    expected["Filename"] = filepath

    # The template options should always be valid
    assert _read_makefile(filepath) == expected


########################################################################################
# Tests for genome specifications


def _write_genome_yaml(tmp_path: Path, genomes: dict[str, str]) -> Path:
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES[:-2]:
            handle.write(read_template(filename))

        handle.write("\nGenomes:\n")
        for key, value in genomes.items():
            handle.write(f"  {key}:\n")
            handle.write(f"    Path: {value}\n")

    return filepath


def test_genomes__no_genomes(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {})
    makefile = _read_makefile(filepath)

    assert makefile["Genomes"] == {}


def test_genomes__simple_genome(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {"Test": "/path/to/genome.fasta"})
    makefile = _read_makefile(filepath)

    assert makefile["Genomes"] == {
        "Test": {"IndexFormat": ".bai", "Name": "Test", "Path": "/path/to/genome.fasta"}
    }


def test_genomes__bad_wildcard_in_name(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {"T*est": "/path/to/genome.fasta"})

    with pytest.raises(MakefileError, match="'*' is not allowed in Genome names"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    (tmp_path / "foo.fasta").write_text(">foo\nACGT\n")
    (tmp_path / "bar.fasta").write_text(">bar\nACGT\n")

    makefile = _read_makefile(filepath)
    assert makefile["Genomes"] == {
        "foo": {"IndexFormat": ".bai", "Name": "foo", "Path": f"{tmp_path}/foo.fasta"},
        "bar": {"IndexFormat": ".bai", "Name": "bar", "Path": f"{tmp_path}/bar.fasta"},
    }


def test_genomes__wildcard_in_name__duplicate_name_1(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(
        tmp_path,
        {
            "foo": "/path/to/genome.fasta",
            "Test*": f"{tmp_path}/*.fasta",
        },
    )

    (tmp_path / "foo.fasta").touch()

    with pytest.raises(MakefileError, match="Multiple genomes named foo"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name__duplicate_name_2(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(
        tmp_path,
        {
            "Test*": f"{tmp_path}/*.fasta",
            "foo": "/path/to/genome.fasta",
        },
    )

    (tmp_path / "foo.fasta").touch()

    with pytest.raises(MakefileError, match="Multiple genomes named foo"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name__no_genomes(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    with pytest.raises(MakefileError, match="Did not find any genomes"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name__invalid_file(tmp_path: Path) -> None:
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    (tmp_path / "foo.fasta").touch()

    with pytest.raises(MakefileError, match="Error reading/indexing FASTA"):
        _read_makefile(filepath)


########################################################################################
# Tests for sample specifications


def _sample_template(filepath: Path) -> dict[str, Any]:
    data = copy.deepcopy(TEMPLATE_BAM)
    data["Filename"] = filepath
    options = data["Options"]

    data["Samples"] = {
        "Sam1": {
            "Lib1": {
                "Run1": [
                    {
                        "Type": "Untrimmed",
                        "Path": ("/path/to/run1.fq.gz", None),
                        "Shortname": "run1.fq.gz",
                        "Options": copy.deepcopy(options),
                    }
                ],
                "Run2": [
                    {
                        "Type": "Untrimmed",
                        "Path": ("/path/to/run2_1.fq.gz", "/path/to/run2_2.fq.gz"),
                        "Shortname": "run2_x.fq.gz",
                        "Options": copy.deepcopy(options),
                    }
                ],
            },
            "Lib2": {
                "Run3": [
                    {
                        "Type": "Untrimmed",
                        "Path": ("/path/to/run3.fq.gz", None),
                        "Shortname": "run3.fq.gz",
                        "Options": copy.deepcopy(options),
                    }
                ]
            },
        }
    }

    return data


def _write_sample_yaml(tmp_path: Path, text: str) -> Path:
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES[:-1]:
            handle.write(read_template(filename))

        handle.write(textwrap.dedent(text))

    return filepath


def test_makefile__sample__basic(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3: /path/to/run3.fq.gz
    """,
    )

    assert _read_makefile(filepath) == _sample_template(filepath)


def test_makefile__sample__group_options(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3: /path/to/run3.fq.gz
          Options:
            AdapterRemoval:
              --adapter1: CATCATDOGCAT
    """,
    )

    expected = _sample_template(filepath)
    for lib, run in [("Lib1", "Run1"), ("Lib1", "Run2"), ("Lib2", "Run3")]:
        options = expected["Samples"]["Sam1"][lib][run][0]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__sample_options(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3: /path/to/run3.fq.gz
            Options:
              AdapterRemoval:
                --adapter1: CATCATDOGCAT
    """,
    )

    expected = _sample_template(filepath)
    for lib, run in [("Lib1", "Run1"), ("Lib1", "Run2"), ("Lib2", "Run3")]:
        options = expected["Samples"]["Sam1"][lib][run][0]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__library_options(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
              Options:
                AdapterRemoval:
                  --adapter1: CATCATDOGCAT
            Lib2:
              Run3: /path/to/run3.fq.gz
    """,
    )

    expected = _sample_template(filepath)
    for lane in ("Run1", "Run2"):
        for record in expected["Samples"]["Sam1"]["Lib1"][lane]:
            record["Options"]["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__lane_options(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2:
                Untrimmed: /path/to/run2_{Pair}.fq.gz
                Options:
                  AdapterRemoval:
                    --adapter1: CATCATDOGCAT
            Lib2:
              Run3: /path/to/run3.fq.gz
    """,
    )

    expected = _sample_template(filepath)
    for record in expected["Samples"]["Sam1"]["Lib1"]["Run2"]:
        record["Options"]["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__mixed_options(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3:
                Untrimmed: /path/to/run3.fq.gz
                Options:
                  AdapterRemoval:
                    --adapter1: DOGCATDOGDOG
            Options:
              AdapterRemoval:
                --adapter1: CATCATDOGCAT
    """,
    )

    expected = _sample_template(filepath)

    for lane in ("Run1", "Run2"):
        for record in expected["Samples"]["Sam1"]["Lib1"][lane]:
            record["Options"]["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    for record in expected["Samples"]["Sam1"]["Lib2"]["Run3"]:
        record["Options"]["AdapterRemoval"]["--adapter1"] = "DOGCATDOGDOG"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__multiple_samples(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
        MoreSamples:
          Sam2:
            Lib2:
              Run3: /path/to/run3.fq.gz
    """,
    )

    expected = _sample_template(filepath)
    expected["Samples"] = {
        "Sam1": {"Lib1": expected["Samples"]["Sam1"]["Lib1"]},
        "Sam2": {"Lib2": expected["Samples"]["Sam1"]["Lib2"]},
    }

    assert _read_makefile(filepath) == expected


def test_makefile__sample__duplicate_sample_names_1(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/trimmed_se.fq.gz
        Others:
          Sam1:
            Lib2:
              Run3: /path/to/run3.fq.gz
    """,
    )

    with pytest.raises(MakefileError, match="Duplicate sample names found"):
        _read_makefile(filepath)


########################################################################################
# Limited options -- can only be modified at certain levels


FEATURES = ("Coverage", "Depths", "mapDamage", "PCRDuplicates", "Summary")


@pytest.mark.parametrize("key", FEATURES)
def test_makefile__sample__features__group(tmp_path: Path, key: str) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
          Options:
            Features:
              {key}: no
    """,
    )

    result = _read_makefile(filepath)
    for record in result["Samples"]["Sam1"]["Lib1"]["Run1"]:
        assert not record["Options"]["Features"][key]


@pytest.mark.parametrize("key", FEATURES)
def test_makefile__sample__features__sample(tmp_path: Path, key: str) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
          Options:
            Features:
              {key}: no
    """,
    )

    result = _read_makefile(filepath)
    for record in result["Samples"]["Sam1"]["Lib1"]["Run1"]:
        assert not record["Options"]["Features"][key]


@pytest.mark.parametrize("key", ["mapDamage", "PCRDuplicates"])
def test_makefile__sample__features__library__allowed(tmp_path: Path, key: str) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Options:
                Features:
                  {key}: no
    """,
    )

    result = _read_makefile(filepath)
    for record in result["Samples"]["Sam1"]["Lib1"]["Run1"]:
        assert not record["Options"]["Features"][key]


@pytest.mark.parametrize("key", ["Coverage", "Depths", "Summary"])
def test_makefile__sample__features__library__prohibited(
    tmp_path: Path, key: str
) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1: /path/to/run1.fq.gz
              Options:
                Features:
                  {key}: yes
    """,
    )

    with pytest.raises(MakefileError, match=f"Cannot override {key!r} at"):
        _read_makefile(filepath)


@pytest.mark.parametrize("key", FEATURES)
def test_makefile__sample__features__lane(tmp_path: Path, key: str) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1:
                Untrimmed: /path/to/run1.fq.gz
                Options:
                  Features:
                    {key}: yes
    """,
    )

    with pytest.raises(MakefileError, match=f"Cannot override {key!r} at"):
        _read_makefile(filepath)


def test_makefile__sample__mapdamage_at_lane_level(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1:
                Untrimmed: /path/to/run1.fq.gz
                Options:
                  mapDamage:
                    --ymax: 17.0
    """,
    )

    with pytest.raises(
        MakefileError,
        match="Cannot set mapDamage options at",
    ):
        _read_makefile(filepath)


########################################################################################
# Pre-trimmed reads


def test_makefile__sample__pretrimmed_reads(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1:
                Single: /path/to/trimmed_se.fq.gz
                Paired: /path/to/trimmed_pe_{Pair}.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3:
                Untrimmed: /path/to/run3.fq.gz
    """,
    )

    expected = _sample_template(filepath)
    expected["Samples"]["Sam1"]["Lib1"]["Run1"] = [
        {
            "Path": ("/path/to/trimmed_se.fq.gz", None),
            "Shortname": "trimmed_se.fq.gz",
            "Type": "Single",
            "Options": expected["Options"],
        },
        {
            "Path": ("/path/to/trimmed_pe_1.fq.gz", "/path/to/trimmed_pe_2.fq.gz"),
            "Shortname": "trimmed_pe_x.fq.gz",
            "Type": "Paired",
            "Options": expected["Options"],
        },
    ]

    assert _read_makefile(filepath) == expected


def test_makefile__sample__pretrimmed_reads__including_untrimmed(
    tmp_path: Path,
) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1:
                Untrimmed: /path/to/trimmed_se.fq.gz
                Paired: /path/to/trimmed_pe_{Pair}.fq.gz
              Run2: /path/to/run2_{Pair}.fq.gz
            Lib2:
              Run3:
                Untrimmed: /path/to/run3.fq.gz
    """,
    )

    with pytest.raises(MakefileError, match="both untrimmed and trimmed reads at"):
        _read_makefile(filepath)


########################################################################################
# Post-processing validation


@pytest.mark.parametrize("offset", [64, "Solexa"])
def test_makefile__validation__trimmed_qualities_must_be_33(
    tmp_path: Path, offset: str | int
) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1:
                Paired: /path/to/trimmed_pe_{{Pair}}.fq.gz
                Options:
                  QualityOffset: {offset}
    """,
    )

    with pytest.raises(MakefileError, match="trimmed data must have quality offset 33"):
        _read_makefile(filepath)


def test_makefile__validation__paired_must_have_key(tmp_path: Path) -> None:
    filepath = _write_sample_yaml(
        tmp_path,
        """
        Samples:
          Sam1:
            Lib1:
              Run1:
                Paired: /path/to/trimmed_pe_Pair.fq.gz
    """,
    )

    with pytest.raises(MakefileError, match="Expected value: a path with a {pair} key"):
        _read_makefile(filepath)
