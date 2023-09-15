import logging
import string

import paleomix.common.yaml as yaml
from paleomix.common.bamfiles import BAM_PLATFORMS
from paleomix.common.makefile import (
    REQUIRED_VALUE,
    And,
    FASTQPath,
    IsBoolean,
    IsFloat,
    IsInt,
    IsListOf,
    IsNone,
    IsStr,
    MakefileError,
    Not,
    Or,
    SpecTree,
    StringStartsWith,
    ValueIn,
    ValuesIntersect,
    ValuesSubsetOf,
    process_makefile,
)


def load_project(filename):
    try:
        with open(filename) as handle:
            data = yaml.safe_load(handle)
    except yaml.YAMLError as error:
        raise MakefileError(error)

    # Fill out any values using user-specified constants; this must be done prior
    # to processing and validating the structure, since the constant names themselves
    # may not be valid values where they are used.
    data = _replace_constants(data)

    data = process_makefile(data, _VALIDATION)

    return _post_process_project(data)


# Valid names for genomes, library, and lanes
_VALID_NAME = And(
    IsStr(),
    ValuesSubsetOf(
        string.ascii_letters + string.digits + "._-",
        description="characters a-z, A-Z, 0-9, '.', '_', and '-' allowed",
    ),
)

# Valid names for samples. Must not overlap with special folders
_VALID_SAMPLE_NAME = And(_VALID_NAME, Not(ValueIn(["reports"])))

# Valid paths for genomes; this is a conservative set of characters chosen to ensure
# that they won't cause problems with poorly written tools and scripts, such as Bash
# scripts that do not quote filenames.
_VALID_GENOME_PATH = And(
    IsStr(),
    Not(ValuesIntersect('\\:?"<>|() \t\n\v\f\r')),
    default=REQUIRED_VALUE,
)

# Values that may used for generic command-line options
_COMMAND_LINE_VALUE = Or(
    IsListOf(IsStr, IsInt, IsFloat), Or(IsStr, IsInt, IsFloat, IsNone)
)

# Command-line options starting with double-dashes; this is done to ensure that only
# one variaint of a command-line option is used (i.e. not both -t and --threads), and
# to make project files more readable.
_LONG_COMMAND_LINE_OPTIONS: SpecTree = {
    StringStartsWith("--"): _COMMAND_LINE_VALUE,
}

# Command-line options starting with single-dashes; should only be used if the tool
# does not support --long-form command-line options.
_SHORT_COMMAND_LINE_OPTIONS: SpecTree = {
    StringStartsWith("-"): _COMMAND_LINE_VALUE,
}


_VALIDATION: SpecTree = {
    "Samples": {
        _VALID_SAMPLE_NAME: {  # Sample
            _VALID_NAME: {  # Library
                _VALID_NAME: FASTQPath(paired_end=True),  # Lane
            },
        },
    },
    "Genome": {
        "Name": _VALID_NAME,
        "Path": _VALID_GENOME_PATH,
    },
    "Settings": {
        "Metadata": {
            "Platform": ValueIn(BAM_PLATFORMS, default=REQUIRED_VALUE),
        },
        #
        "Constants": {
            IsStr: _COMMAND_LINE_VALUE,
        },
        # Common options passe to the JRE
        "JavaOptions": {
            StringStartsWith("-X"): IsNone,
        },
        "Preprocessing": {
            # Quality metrics performed prior to pre-processsing
            "FastQC": _LONG_COMMAND_LINE_OPTIONS,
        },
        "ReadMapping": {
            # FIXME: Should be grouped with other mapping-step programs (fixmate, etc.)
            "BWAMem": _SHORT_COMMAND_LINE_OPTIONS,
            "PCRDuplicates": {
                "mode": ValueIn(("mark", "filter", "skip"), default="mark"),
            },
            "BaseRecalibrator": {
                StringStartsWith("--"): _COMMAND_LINE_VALUE,
                "--known-sites": IsStr(default=REQUIRED_VALUE),
            },
            "ApplyBQSR": _LONG_COMMAND_LINE_OPTIONS,
        },
        "Genotyping": {
            "HaplotypeCaller": _LONG_COMMAND_LINE_OPTIONS,
            "GenotypeGVCFs": _LONG_COMMAND_LINE_OPTIONS,
            "VariantRecalibrator": {
                "Enabled": IsBoolean(default=True),
                "INDEL": _LONG_COMMAND_LINE_OPTIONS,
                "SNP": _LONG_COMMAND_LINE_OPTIONS,
            },
            "ApplyVQSR": {
                "Enabled": IsBoolean(default=True),
                "INDEL": _LONG_COMMAND_LINE_OPTIONS,
                "SNP": _LONG_COMMAND_LINE_OPTIONS,
            },
        },
    },
}


########################################################################################
# Pre- and post-processing of projects


def _replace_constants(data):
    settings = data.get("Settings", {})
    constants = settings.pop("Constants", {})
    if not isinstance(constants, dict):
        # Invalid constants will be caught during global processing
        return

    # Validate user-supplied constants; the structure is replicated to ensure that
    # error messages give the correct path into the settings structure
    process_makefile(
        {"Settings": {"Constants": constants}},
        {"Settings": {"Constants": {IsStr: _COMMAND_LINE_VALUE}}},
    )
    constants = {("${%s}" % (key)): value for key, value in constants.items()}
    if not constants:
        return

    data["Settings"] = _replace_constants_recursive(settings, constants)

    return data


def _replace_constants_recursive(data, constants):
    if isinstance(data, dict):
        return {
            _replace_constants_recursive(key, constants): _replace_constants_recursive(
                value, constants
            )
            for key, value in data.items()
        }
    elif isinstance(data, list):
        return [_replace_constants_recursive(value, constants) for value in data]
    elif isinstance(data, str):
        return constants.get(data, data)
    else:
        return data


def _post_process_project(data):
    post_processing_steps = [
        _process_fastq_paths,
        _process_gatk_resources,
        _process_recalibrator_settings,
    ]

    any_errors = False
    for mangle_func in post_processing_steps:
        if not mangle_func(data):
            any_errors = True

    if any_errors:
        raise MakefileError("Invalid settings in project")

    return data


def _process_fastq_paths(data):
    for libraries in data["Samples"].values():
        for lanes in libraries.values():
            for lane, filename in lanes.items():
                lanes[lane] = {
                    1: FASTQPath.format(filename, 1),
                    2: FASTQPath.format(filename, 2),
                }

    return True


def _process_gatk_resources(data):
    any_errors = False
    log = logging.getLogger(__name__)
    # Validate resource files used by VariantRecalibrator
    for mode, options in data["Settings"]["Genotyping"]["VariantRecalibrator"].items():
        if mode not in ("SNP", "INDEL"):
            continue

        found_any = False
        for key, value in options.items():
            # Keys take the form `--resource` and `--resource:${settings}`
            if key == "--resource" or key.startswith("--resource:"):
                if not (isinstance(value, str) and value.endswith(".vcf.gz")):
                    any_errors = True
                    log.error(
                        "Invalid %s path for %s recalibrator; expected a bgzip "
                        "compressed VCF file (*.vcf.gz) but found %r",
                        key,
                        mode,
                        value,
                    )

                found_any = True

        if not found_any:
            log.error("No --resource files specified for %s recalibrator", mode)
            any_errors = True

    return not any_errors


def _process_recalibrator_settings(data):
    # ApplyVQSR relies on files created by the VariantRecalibrator step. For simplicity
    # we require that both be enabled for ApplyVQSR to be carried out. It would be
    # possible to run ApplyVQSR with the prior step disabled, provided that the models
    # had been built, but this would result in unexpected behavior if upstream files
    # changed (models not beeing updated).
    genotyping = data["Settings"]["Genotyping"]
    if not genotyping["VariantRecalibrator"]["Enabled"]:
        if genotyping["ApplyVQSR"]["Enabled"]:
            genotyping["ApplyVQSR"]["Enabled"] = False
            log = logging.getLogger(__name__)
            log.warning(
                "GATK variant recalibration (ApplyVQSR) is enabled, but model building "
                "(VariantRecalibrator) is disabled. Variant recalibration will NOT be "
                "performed!"
            )

    return True
