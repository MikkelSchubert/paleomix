#!/usr/bin/python3
import os

from paleomix.node import CommandNode
from paleomix.common.fileutils import describe_files
from paleomix.tools import factory


class DetectInputDuplicationNode(CommandNode):
    def __init__(self, input_files, output_file, dependencies=()):
        command = factory.new("dupcheck")
        command.set_option("--pipeline-output", os.path.abspath(output_file))
        command.add_multiple_values(input_files)
        command.set_kwargs(OUT_STDOUT=output_file)

        CommandNode.__init__(
            self,
            description="detecting duplicate input in %s"
            % (describe_files(input_files)),
            command=command.finalize(),
            dependencies=dependencies,
        )


class ValidateFASTQFilesNode(CommandNode):
    def __init__(
        self, input_files, output_file, offset, collapsed=False, dependencies=()
    ):
        command = factory.new(":validate_fastq")
        command.set_option("--offset", offset)
        if collapsed:
            command.set_option("--collapsed")
        command.add_multiple_values(input_files)
        command.set_kwargs(OUT_STDOUT=output_file)

        CommandNode.__init__(
            self,
            description="validating %s" % (describe_files(input_files),),
            command=command.finalize(),
            dependencies=dependencies,
        )


class ValidateFASTAFilesNode(CommandNode):
    def __init__(self, input_file, output_file, dependencies=()):
        command = factory.new(":validate_fasta")
        command.add_value("%(IN_FASTA)s")
        command.set_kwargs(IN_FASTA=input_file, OUT_STDOUT=output_file)

        CommandNode.__init__(
            self,
            description="validating %s" % (input_file,),
            command=command.finalize(),
            dependencies=dependencies,
        )
