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
from paleomix.nodes.picard import ValidateBAMNode
from paleomix.nodes.samtools import BAMIndexNode


def index_and_validate_bam(
    config,
    prefix,
    node,
    log_file=None,
    create_index=True,
    validation_levels=("full",),
):
    input_file, index_file = _get_input_files(node, prefix["IndexFormat"])
    if not index_file and create_index:
        node = BAMIndexNode(
            infile=input_file, index_format=prefix["IndexFormat"], dependencies=node
        )
        (index_file,) = node.output_files

    # Only enable validation if the user picked a corresponding level
    if config.validation not in validation_levels:
        return node

    ignored_checks = [
        # Ignored since we may filter out misses and low-quality hits during
        # mapping, which leads to a large proportion of missing PE mates.
        "MATE_NOT_FOUND",
        # Ignored due to high rate of false positives for lanes with few hits,
        # where high-quality reads may cause mis-identification of qualities
        "INVALID_QUALITY_FORMAT",
    ]

    return ValidateBAMNode(
        config=config,
        input_bam=input_file,
        input_index=index_file,
        ignored_checks=ignored_checks,
        big_genome_mode=prefix["IndexFormat"] == ".csi",
        output_log=log_file,
        dependencies=node,
    )


def _get_input_files(node, index_format):
    index_filename = None
    input_filename = None
    for filename in node.output_files:
        if filename.lower().endswith(index_format):
            index_filename = True
        elif filename.lower().endswith(".bam"):
            input_filename = filename

    return input_filename, index_filename
