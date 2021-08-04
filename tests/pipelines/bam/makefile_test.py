import copy
import textwrap

import pytest

from paleomix.common.makefile import MakefileError
from paleomix.pipelines.bam.makefile import read_makefiles
from paleomix.resources import template


def _read_makefile(filepath):
    # The template options should always be valid
    (result,) = read_makefiles([filepath])

    # TODO: Remove once no longer needed for backwards compatibility
    del result["Targets"]

    return result


########################################################################################
# Test the validity of the makefile templates

# Options specified in a basic YAML file
TEMPLATE_BAM_OPTIONS = {
    "AdapterRemoval": {
        "--collapse": None,
        "--minlength": 25,
        "--trimns": None,
        "--trimqualities": None,
    },
    "Aligners": {
        "Program": "BWA",
        "BWA": {
            "Algorithm": "backtrack",
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
    "Platform": "Illumina",
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
                "NAME_OF_LANE::001": {
                    "Data": {"SE": "PATH_WITH_WILDCARDS"},
                    "Options": copy.deepcopy(TEMPLATE_BAM_OPTIONS),
                    "Tags": {
                        "DS": "PATH_WITH_WILDCARDS",
                        "Folder": "PATH_WITH_WILDCARDS",
                        "ID": "NAME_OF_LIBRARY",
                        "LB": "NAME_OF_LIBRARY",
                        "PG": "BWA",
                        "PL": "ILLUMINA",
                        "PU": "NAME_OF_LANE",
                        "SM": "NAME_OF_SAMPLE",
                        "Target": "NAME_OF_SAMPLE",
                        "Group": "Samples",
                    },
                    "Type": "Untrimmed",
                }
            }
        }
    },
}


def _trim_template():
    data = copy.deepcopy(TEMPLATE_BAM)
    options = [
        data,
        data["Samples"]["NAME_OF_SAMPLE"]["NAME_OF_LIBRARY"]["NAME_OF_LANE::001"],
    ]

    del data["Genomes"]["NAME_OF_PREFIX"]

    for dd in options:
        del dd["Options"]["Aligners"]["Bowtie2"]["--very-sensitive"]
        del dd["Options"]["mapDamage"]["--downsample"]

    return data


TEMPLATE_TRIM = _trim_template()


def test_basic__read_makefiles__trim_template(tmp_path):
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        # Minimal template containing options required for trim pipeline
        handle.write(template("bam_head.yaml"))
        # Example samples (not used when auto generating a YAML file)
        handle.write(template("bam_samples.yaml"))

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


def test_basic__read_makefiles__bam_template(tmp_path):
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES:
            handle.write(template(filename))

    expected = copy.deepcopy(TEMPLATE_BAM)
    expected["Filename"] = filepath

    # The template options should always be valid
    assert _read_makefile(filepath) == expected


########################################################################################
# Tests for genome specifications


def _genome_template():
    template = copy.deepcopy(TEMPLATE_BAM)
    template.pop("Samples")
    template["Genomes"].pop("NAME_OF_PREFIX")

    return template


GENOMES_TEMPLATE = _genome_template()


def _write_genome_yaml(tmp_path, genomes):
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES[:-2]:
            handle.write(template(filename))

        handle.write("\nGenomes:\n")
        for key, value in genomes.items():
            handle.write(f"  {key}:\n")
            handle.write(f"    Path: {value}\n")

    return filepath


def test_genomes__no_genomes(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {})
    makefile = _read_makefile(filepath)

    assert makefile["Genomes"] == {}


def test_genomes__simple_genome(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {"Test": "/path/to/genome.fasta"})
    makefile = _read_makefile(filepath)

    assert makefile["Genomes"] == {
        "Test": {"IndexFormat": ".bai", "Name": "Test", "Path": "/path/to/genome.fasta"}
    }


def test_genomes__bad_wildcard_in_name(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {"T*est": "/path/to/genome.fasta"})

    with pytest.raises(MakefileError, match="'*' is not allowed in Genome names"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    (tmp_path / "foo.fasta").write_text(">foo\nACGT\n")
    (tmp_path / "bar.fasta").write_text(">bar\nACGT\n")

    makefile = _read_makefile(filepath)
    assert makefile["Genomes"] == {
        "foo": {"IndexFormat": ".bai", "Name": "foo", "Path": f"{tmp_path}/foo.fasta"},
        "bar": {"IndexFormat": ".bai", "Name": "bar", "Path": f"{tmp_path}/bar.fasta"},
    }


def test_genomes__wildcard_in_name__duplicate_name_1(tmp_path):
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


def test_genomes__wildcard_in_name__duplicate_name_2(tmp_path):
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


def test_genomes__wildcard_in_name__no_genomes(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    with pytest.raises(MakefileError, match="Did not find any genomes"):
        _read_makefile(filepath)


def test_genomes__wildcard_in_name__invalid_file(tmp_path):
    filepath = _write_genome_yaml(tmp_path, {"Test*": f"{tmp_path}/*.fasta"})

    (tmp_path / "foo.fasta").touch()

    with pytest.raises(MakefileError, match="Error reading/indexing FASTA"):
        _read_makefile(filepath)


########################################################################################
# Tests for sample specifications


def _sample_template(filepath):
    data = copy.deepcopy(TEMPLATE_BAM)
    data["Filename"] = filepath
    options = data["Options"]

    data["Samples"] = {
        "Sam1": {
            "Lib1": {
                "Run1::001": {
                    "Type": "Untrimmed",
                    "Data": {"SE": "/path/to/run1.fq.gz"},
                    "Options": copy.deepcopy(options),
                    "Tags": {
                        "Target": "Sam1",
                        "Group": "Samples",
                        "ID": "Lib1",
                        "SM": "Sam1",
                        "LB": "Lib1",
                        "PU": "Run1",
                        "DS": "/path/to/run1.fq.gz",
                        "Folder": "run1.fq.gz",
                        "PG": "BWA",
                        "PL": "ILLUMINA",
                    },
                },
                "Run2::001": {
                    "Type": "Untrimmed",
                    "Data": {
                        "PE_1": "/path/to/run2_1.fq.gz",
                        "PE_2": "/path/to/run2_2.fq.gz",
                    },
                    "Options": copy.deepcopy(options),
                    "Tags": {
                        "Target": "Sam1",
                        "Group": "Samples",
                        "ID": "Lib1",
                        "SM": "Sam1",
                        "LB": "Lib1",
                        "PU": "Run2",
                        "DS": "/path/to/run2_[12].fq.gz",
                        "Folder": "run2_x.fq.gz",
                        "PG": "BWA",
                        "PL": "ILLUMINA",
                    },
                },
            },
            "Lib2": {
                "Run3::001": {
                    "Type": "Untrimmed",
                    "Data": {"SE": "/path/to/run3.fq.gz"},
                    "Options": copy.deepcopy(options),
                    "Tags": {
                        "Target": "Sam1",
                        "Group": "Samples",
                        "ID": "Lib2",
                        "SM": "Sam1",
                        "LB": "Lib2",
                        "PU": "Run3",
                        "DS": "/path/to/run3.fq.gz",
                        "Folder": "run3.fq.gz",
                        "PG": "BWA",
                        "PL": "ILLUMINA",
                    },
                }
            },
        }
    }

    return data


def _write_sample_yaml(tmp_path, text):
    filepath = tmp_path / "temp.yaml"
    with filepath.open("wt") as handle:
        for filename in TEMPLATE_FILES[:-1]:
            handle.write(template(filename))

        handle.write(textwrap.dedent(text))

    return filepath


def test_makefile__sample__basic(tmp_path):
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

    expected = _sample_template(filepath)

    assert _read_makefile(filepath) == expected


def test_makefile__sample__group_options(tmp_path):
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
    for (lib, run) in [("Lib1", "Run1"), ("Lib1", "Run2"), ("Lib2", "Run3")]:
        options = expected["Samples"]["Sam1"][lib][f"{run}::001"]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__sample_options(tmp_path):
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
    for (lib, run) in [("Lib1", "Run1"), ("Lib1", "Run2"), ("Lib2", "Run3")]:
        options = expected["Samples"]["Sam1"][lib][f"{run}::001"]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__library_options(tmp_path):
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
    for (lib, run) in [("Lib1", "Run1"), ("Lib1", "Run2")]:
        options = expected["Samples"]["Sam1"][lib][f"{run}::001"]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__lane_options(tmp_path):
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
    options = expected["Samples"]["Sam1"]["Lib1"][f"Run2::001"]["Options"]
    options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__mixed_options(tmp_path):
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
    for (lib, run) in [("Lib1", "Run1"), ("Lib1", "Run2")]:
        options = expected["Samples"]["Sam1"][lib][f"{run}::001"]["Options"]
        options["AdapterRemoval"]["--adapter1"] = "CATCATDOGCAT"

    options = expected["Samples"]["Sam1"]["Lib2"]["Run3::001"]["Options"]
    options["AdapterRemoval"]["--adapter1"] = "DOGCATDOGDOG"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__multiple_samples(tmp_path):
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

    tags = expected["Samples"]["Sam2"]["Lib2"]["Run3::001"]["Tags"]
    tags["Group"] = "MoreSamples"
    tags["SM"] = "Sam2"
    tags["Target"] = "Sam2"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__duplicate_sample_names_1(tmp_path):
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

    with pytest.raises(MakefileError, match="Multiple samples named"):
        _read_makefile(filepath)


########################################################################################
# Limited options -- can only be modified at certain levels


FEATURES = ("Coverage", "Depths", "mapDamage", "PCRDuplicates", "Summary")


@pytest.mark.parametrize("key", FEATURES)
def test_makefile__sample__features__group(tmp_path, key):
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

    result = _read_makefile(filepath)
    options = result["Samples"]["Sam1"]["Lib1"]["Run1::001"]["Options"]
    assert options["Features"][key] == True


@pytest.mark.parametrize("key", FEATURES)
def test_makefile__sample__features__sample(tmp_path, key):
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

    result = _read_makefile(filepath)
    options = result["Samples"]["Sam1"]["Lib1"]["Run1::001"]["Options"]
    assert options["Features"][key] == True


@pytest.mark.parametrize("key", ("mapDamage", "PCRDuplicates"))
def test_makefile__sample__features__library__allowed(tmp_path, key):
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

    result = _read_makefile(filepath)
    options = result["Samples"]["Sam1"]["Lib1"]["Run1::001"]["Options"]
    assert options["Features"][key] == True


@pytest.mark.parametrize("key", ("Coverage", "Depths", "Summary"))
def test_makefile__sample__features__library__prohibited(tmp_path, key):
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
def test_makefile__sample__features__lane(tmp_path, key):
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


def test_makefile__sample__mapdamage_at_lane_level(tmp_path):
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
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
        match=f"Cannot set mapDamage options for individual lanes",
    ):
        _read_makefile(filepath)


########################################################################################
# Pre-trimmed reads


def test_makefile__sample__pretrimmed_reads(tmp_path):
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
    # No auto-renumbering of split lanes
    run = expected["Samples"]["Sam1"]["Lib1"].pop("Run1::001")
    expected["Samples"]["Sam1"]["Lib1"]["Run1"] = run
    # Minor differences in structure
    run["Type"] = "Trimmed"
    run["Data"] = {
        "Single": "/path/to/trimmed_se.fq.gz",
        "Paired": "/path/to/trimmed_pe_{Pair}.fq.gz",
    }
    run["Tags"]["DS"] = "NA"
    run["Tags"]["Folder"] = "NA"

    assert _read_makefile(filepath) == expected


def test_makefile__sample__pretrimmed_reads__including_untrimmed(tmp_path):
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


@pytest.mark.parametrize("offset", (64, "Solexa"))
def test_makefile__validation__trimmed_qualities_must_be_33(tmp_path, offset):
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


def test_makefile__validation__paired_must_have_key(tmp_path):
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1:
                Paired: /path/to/trimmed_pe_Pair.fq.gz
    """,
    )

    with pytest.raises(MakefileError, match="Paired pre-trimmed reads specified, but"):
        _read_makefile(filepath)


def test_makefile__validation__single_must_not_have_key(tmp_path):
    filepath = _write_sample_yaml(
        tmp_path,
        f"""
        Samples:
          Sam1:
            Lib1:
              Run1:
                Single: /path/to/trimmed_pe_{{Pair}}.fq.gz
    """,
    )

    with pytest.raises(MakefileError, match="read-type is not Paired"):
        _read_makefile(filepath)
