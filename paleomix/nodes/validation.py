#!/usr/bin/python3
from paleomix.atomiccmd.command2 import InputFile
from paleomix.common.fileutils import describe_files
from paleomix.node import CommandNode
from paleomix.tools import factory


class ValidateFASTQFilesNode(CommandNode):
    def __init__(
        self, input_files, output_file, offset, collapsed=False, dependencies=()
    ):
        command = factory.new(
            [":validate_fastq", "--offset", offset],
            stdout=output_file,
        )

        if collapsed:
            command.append("--collapsed")

        for filename in input_files:
            command.append(InputFile(filename))

        CommandNode.__init__(
            self,
            description="validating %s" % (describe_files(input_files),),
            command=command,
            dependencies=dependencies,
        )


class ValidateFASTAFilesNode(CommandNode):
    def __init__(self, input_file, output_file, dependencies=()):
        command = factory.new(
            [":validate_fasta", InputFile(input_file)],
            stdout=output_file,
        )

        CommandNode.__init__(
            self,
            description="validating %s" % (input_file,),
            command=command,
            dependencies=dependencies,
        )
