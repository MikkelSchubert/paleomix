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
import os

from paleomix.common.utilities import safe_coerce_to_tuple
from paleomix.nodes.commands import FilterCollapsedBAMNode
from paleomix.nodes.mapdamage import (
    MapDamageModelNode,
    MapDamagePlotNode,
    MapDamageRescaleNode,
)
from paleomix.nodes.picard import MarkDuplicatesNode
from paleomix.nodes.validation import DetectInputDuplicationNode
from paleomix.pipelines.bam.nodes import index_and_validate_bam


class Library:
    """Represents a single library in a BAM pipeline.

    Is reponsible for aggregating per-lane BAMS, removal of PCR duplicates,
    rescaling of quality-scores using mapDamage, as well as running mapDamage
    for QC purposes.

    Properties:
      name      -- Name of the libray (as specified in makefile)
      lanes     -- Tuple of lanes assosisated with the library
      options   -- Makefile options that apply to the current library
      folder    -- Folder containing files assosisated with library. Is used as
                   a prefix for files generated by this class.
      bams      -- Dictionary of BAM filenames -> nodes, for each BAM generated by
                   the Library class. Depending on options, this may either be newly
                   generated files, or the files produced by Lanes.
    """

    def __init__(self, config, target, prefix, lanes, name):
        self.name = name
        self.lanes = safe_coerce_to_tuple(lanes)
        self.options = lanes[0].options
        self.folder = os.path.dirname(os.path.dirname(self.lanes[0].folder))

        assert all(
            (self.folder == os.path.dirname(os.path.dirname(lane.folder)))
            for lane in self.lanes
        )
        assert all((self.options == lane.options) for lane in self.lanes)

        lane_bams = self._collect_bams_by_type(self.lanes)

        pcr_duplicates = self.options["Features"]["PCRDuplicates"]
        if pcr_duplicates:
            # pcr_duplicates may be "mark" or any trueish value
            lane_bams = self._remove_pcr_duplicates(
                config, prefix, lane_bams, pcr_duplicates
            )

        # At this point we no longer need to differentiate between read types
        files_and_nodes = self._collect_files_and_nodes(lane_bams)

        # Collect output bams, possible following rescaling
        self.bams, mapdamage_nodes = self._build_mapdamage_nodes(
            config, target, prefix, files_and_nodes
        )

        nodes = [self._build_dataduplication_node(lane_bams)]
        nodes.extend(mapdamage_nodes)

        self.nodes = tuple(nodes)

    @classmethod
    def _collect_bams_by_type(cls, lanes):
        bams = {}
        for lane in lanes:
            for key, files in lane.bams.items():
                key = "collapsed" if (key == "Collapsed") else "normal"
                bams.setdefault(key, {}).update(files)

        return bams

    @classmethod
    def _collect_files_and_nodes(cls, bams):
        files_and_nodes = {}
        for dd in bams.values():
            files_and_nodes.update(dd)
        return files_and_nodes

    def _remove_pcr_duplicates(self, config, prefix, bams, strategy):
        rmdup_cls = {"collapsed": FilterCollapsedBAMNode, "normal": MarkDuplicatesNode}

        keep_duplicates = False
        if isinstance(strategy, str) and (strategy.lower() == "mark"):
            keep_duplicates = True

        results = {}
        for key, files_and_nodes in bams.items():
            output_filename = self.folder + ".rmdup.%s.bam" % key
            node = rmdup_cls[key](
                config=config,
                input_bams=list(files_and_nodes.keys()),
                output_bam=output_filename,
                keep_dupes=keep_duplicates,
                dependencies=list(files_and_nodes.values()),
            )

            # Indexing is required if we wish to calulate per-region statistics
            validated_node = index_and_validate_bam(
                config=config,
                prefix=prefix,
                node=node,
                create_index=bool(prefix.get("RegionsOfInterest")),
            )

            results[key] = {output_filename: validated_node}
        return results

    def _build_mapdamage_nodes(self, config, target, prefix, files_and_nodes):
        # Messing with these does not cause the pipeline to re-do other stuff
        destination = os.path.join(
            config.destination, "%s.%s.mapDamage" % (target, prefix["Name"]), self.name
        )

        run_type = self.options["Features"]["mapDamage"]
        if run_type == "rescale":
            return self._mapdamage_rescale(
                config=config,
                destination=destination,
                prefix=prefix,
                files_and_nodes=files_and_nodes,
            )

        elif run_type == "model":
            # Run of mapDamage including both plots and damage models
            node = self._mapdamage_model(
                destination=destination,
                prefix=prefix,
                files_and_nodes=files_and_nodes,
            )

            return files_and_nodes, (node,)
        elif run_type in ("plot", True):
            # Basic run of mapDamage, only generates plots / tables
            node = self._mapdamage_plot(
                destination=destination,
                prefix=prefix,
                files_and_nodes=files_and_nodes,
            )

            return files_and_nodes, (node,)
        else:
            assert not run_type, run_type
            return files_and_nodes, ()

    def _mapdamage_plot(self, destination, prefix, files_and_nodes):
        title = "mapDamage plot for library %r" % (self.name,)

        return MapDamagePlotNode(
            reference=prefix["Path"],
            input_files=list(files_and_nodes),
            output_directory=destination,
            title=title,
            options=self.options["mapDamage"],
            dependencies=files_and_nodes.values(),
        )

    def _mapdamage_model(self, destination, prefix, files_and_nodes):
        # Generates basic plots / table files
        plot = self._mapdamage_plot(
            destination=destination,
            prefix=prefix,
            files_and_nodes=files_and_nodes,
        )

        # Builds model of post-mortem DNA damage
        return MapDamageModelNode(
            reference=prefix["Reference"],
            directory=destination,
            options=self.options["mapDamage"],
            dependencies=plot,
        )

    def _mapdamage_rescale(self, config, destination, prefix, files_and_nodes):
        model = self._mapdamage_model(
            destination=destination,
            prefix=prefix,
            files_and_nodes=files_and_nodes,
        )

        # Rescales BAM quality scores using model built above
        input_files = list(files_and_nodes)
        output_filename = self.folder + ".rescaled.bam"

        scale = MapDamageRescaleNode(
            reference=prefix["Reference"],
            input_files=input_files,
            output_file=output_filename,
            directory=destination,
            options=self.options["mapDamage"],
            dependencies=model,
        )

        # Grab indexing and validation nodes, required by ROIs
        validate = index_and_validate_bam(
            config=config,
            prefix=prefix,
            node=scale,
            create_index=bool(prefix.get("RegionsOfInterest")),
        )

        return {output_filename: validate}, (model,)

    def _build_dataduplication_node(self, bams):
        files_and_nodes = self._collect_files_and_nodes(bams)
        output_file = self.folder + ".duplications_checked"

        return DetectInputDuplicationNode(
            input_files=list(files_and_nodes),
            output_file=output_file,
            dependencies=list(files_and_nodes.values()),
        )
