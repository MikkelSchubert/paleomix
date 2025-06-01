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
import logging
import os
import shutil
import sys

from paleomix.common.fileutils import copy_file


def _parse_resource(resource):
    *modules, filename = resource.split(os.path.sep)
    return ".".join(modules), filename


if sys.version_info < (3, 7):
    from pkg_resources import cleanup_resources, resource_filename

    class _ResourceContext:
        def __init__(self, resource):
            self._resource = resource
            self._source = None

        def __enter__(self):
            module, filename = _parse_resource(self._resource)
            self._source = resource_filename(module, filename)
            return self

        def __exit__(self, type, value, traceback):
            cleanup_resources()

    def access(resource):
        return _ResourceContext(resource)

elif sys.version_info < (3, 9):
    from importlib import resources

    def access(resource):
        module, filename = _parse_resource(resource)
        return resources.path("paleomix.resources.{}".format(module), filename)

else:
    from importlib import resources

    def _as_traversable(resource):
        return resources.files("paleomix").joinpath("resources").joinpath(resource)

    def access(resource):
        return resources.as_file(_as_traversable(resource))


def copy_resource(resource, destination):
    with access(resource) as path:
        if path.is_dir():
            shutil.copytree(path, destination)
        else:
            copy_file(path, destination)


def add_copy_example_command(subparsers):
    parser = subparsers.add_parser("example", help="Create example project")

    parser.add_argument(
        "destination",
        default=".",
        nargs="?",
        help="Destination folder for example data.",
    )


def copy_example(tool, args):
    """Command-line interface to copy a folder containing example data to a
    folder specified by the user. Arguments are a tool name (e.g.
    'bam_pipeline'), and any command-line options specified by the user;
    returns 0 on success, or 1 on errors.
    """
    log = logging.getLogger(__name__)
    destination = os.path.join(args.destination, tool)
    log.info("Copying example project to %r", args.destination)

    if os.path.exists(destination):
        log.error("Example folder already exists at %r", destination)
        return 1

    copy_resource(os.path.join("examples", tool), destination)

    log.info("Sucessfully saved example in %r", destination)

    return 0
