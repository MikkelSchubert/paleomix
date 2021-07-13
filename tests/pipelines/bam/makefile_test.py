import copy

from paleomix.pipelines.bam.makefile import read_makefiles
from paleomix.resources import template

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
                    },
                    "Type": "Raw",
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
    (result,) = read_makefiles([filepath])

    # TODO: Remove once no longer needed for backwards compatibility
    del result["Targets"]
    del result["Prefixes"]

    assert result == expected


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
    (result,) = read_makefiles([filepath])

    # TODO: Remove once no longer needed for backwards compatibility
    del result["Targets"]
    del result["Prefixes"]

    assert result == expected
