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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#
from nose.tools import assert_equals

from pypeline.common.fileutils import *



################################################################################
################################################################################
## Tests for 'add_postfix'

def test_add_postfix__():
    assert False, "TODO"




################################################################################
################################################################################
## Tests for 'swap_ext'

def test_swap_ext__has_ext_vs_empty_ext():
    assert_equals(swap_ext("name.foo", ""), "name")

def test_swap_ext__empty_ext_vs_empty_ext():
    assert_equals(swap_ext("name", ""), "name")

def test_swap_ext__has_ext_vs_dot_ext():
    assert_equals(swap_ext("name.foo", "."), "name")

def test_swap_ext__dot_ext_vs_dot_ext():
    assert_equals(swap_ext("name.", "."), "name")

# FIXME def test_swap_ext__


################################################################################
################################################################################
## Tests for 'reroot_path'

def test_reroot_path__():
    assert False, "TODO"




################################################################################
################################################################################
## Tests for 'missing_files'

def test_missing_files__file_exists():
    assert_equals(missing_files(["tests/data/empty_file_1"]), [])

def test_missing_files__file_doesnt_exist():
    assert_equals(missing_files(["tests/data/missing_file_1"]), 
                 ["tests/data/missing_file_1"])

def test_missing_files__mixed_files():
    files = ["tests/data/missing_file_1",
             "tests/data/empty_file_1"]
    result = ["tests/data/missing_file_1"]

    assert_equals(missing_files(files), result)




################################################################################
################################################################################
## Tests for 'modified_after'

def test_modified_after__():
    assert False, "TODO"




################################################################################
################################################################################
## Tests for 'is_executable'

def test_is_executable__full_path__is_executable():
    assert is_executable("/bin/ls")

def test_is_executable__full_path__is_non_executable():
    assert not is_executable("/etc/fstab")

def test_is_executable__rel_path__is_executable():
    assert is_executable("tests/unit/run.sh")

def test_is_executable__rel_path__is_non_executable():
    assert not is_executable("tests/data/empty_file_1")




################################################################################
################################################################################
## Tests for 'executable_exists'

def test_executable_exists__executable():
    assert executable_exists("ls")

def test_executable_exists__non_executable():
    assert not executable_exists("lsxxxx")

def test_executable_exists__full_path__is_executable():
    assert executable_exists("/bin/ls")

def test_executable_exists__full_path__is_non_executable():
    assert not executable_exists("/etc/fstab")

def test_executable_exists__rel_path__is_executable():
    assert executable_exists("tests/unit/run.sh")

def test_executable_exists__rel_path__is_non_executable():
    assert not executable_exists("tests/data/empty_file_1")




################################################################################
################################################################################
## Tests for 'missing_executables'

def test_missing_executables__executable():
    assert_equals(missing_executables(["ls"]), [])

def test_missing_executables__non_executable():
    assert_equals(missing_executables(["lsxxxx"]), ["lsxxxx"])

def test_missing_executables__mixed():
    assert_equals(missing_executables(["lsxxxx", "ls"]), ["lsxxxx"])




################################################################################
################################################################################
## Tests for 'move_file'

# FIXME


################################################################################
################################################################################
## Tests for 'copy_file'


# FIXME
