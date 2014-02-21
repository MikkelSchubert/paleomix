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
from pypeline.node import \
    CommandNode, \
    Node
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder
from pypeline.nodes.picard import \
    concatenate_input_bams
from pypeline.common.fileutils import \
    reroot_path, \
    move_file, \
    describe_files
from pypeline.common.utilities import \
    safe_coerce_to_tuple

import pypeline.tools.bam_stats.coverage as coverage


class CoverageNode(CommandNode):
    def __init__(self, config, target_name, input_files, output_file,
                 regions_file=None, dependencies=()):
        kwargs = {"OUT_FILE": output_file}
        input_files = safe_coerce_to_tuple(input_files)

        if len(input_files) > 1:
            if regions_file:
                raise ValueError("Coverage for regions require single, "
                                 "indexed input BAM file.")
            cat_cmds, cat_obj = concatenate_input_bams(config, input_files)
            kwargs["IN_STDIN"] = cat_obj
            input_argument = "-"  # Read from STDIN
        else:
            cat_cmds = []
            kwargs["IN_FILE"] = input_files[0]
            input_argument = "%(IN_FILE)s"

        call = ['bam_coverage', input_argument, '%(OUT_FILE)s',
                '--target-name', target_name]
        builder = AtomicCmdBuilder(call, **kwargs)
        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        command = ParallelCmds(cat_cmds + [builder.finalize()])
        CommandNode.__init__(self,
                             command=command,
                             description="<Coverage: %s -> '%s'>"
                                         % (describe_files(input_files),
                                            output_file),
                             dependencies=dependencies)


class MergeCoverageNode(Node):
    def __init__(self, input_files, output_file, dependencies=()):
        self._output_file = output_file

        Node.__init__(self,
                      description="<MergeCoverage: '%s' -> '%s'>"
                      % (describe_files(input_files), self._output_file),
                      input_files=input_files,
                      output_files=self._output_file,
                      dependencies=dependencies)

    def _run(self, _config, temp):
        table = {}
        for filename in self.input_files:
            coverage.read_table(table, filename)

        coverage.write_table(table, reroot_path(temp, self._output_file))
        move_file(reroot_path(temp, self._output_file), self._output_file)


class DepthHistogramNode(CommandNode):
    def __init__(self, config, target_name, input_files, output_file,
                 regions_file=None, dependencies=()):
        kwargs = {"OUT_FILE": output_file}
        input_files = safe_coerce_to_tuple(input_files)

        if len(input_files) > 1:
            if regions_file:
                raise ValueError("DepthHistogram for regions require single, "
                                 "indexed input BAM file.")
            cat_cmds, cat_obj = concatenate_input_bams(config, input_files)
            kwargs["IN_STDIN"] = cat_obj
            input_argument = "-"  # Read from STDIN
        else:
            cat_cmds = []
            kwargs["IN_FILE"] = input_files[0]
            input_argument = "%(IN_FILE)s"

        call = ['bam_depths', input_argument, '%(OUT_FILE)s',
                '--target-name', target_name]
        builder = AtomicCmdBuilder(call, **kwargs)
        if regions_file:
            builder.set_option('--regions-file', '%(IN_REGIONS)s')
            builder.set_kwargs(IN_REGIONS=regions_file)

        command = ParallelCmds(cat_cmds + [builder.finalize()])
        CommandNode.__init__(self,
                             command=command,
                             description="<DepthHistogram: %s -> '%s'>"
                                         % (describe_files(input_files),
                                            output_file),
                             dependencies=dependencies)
