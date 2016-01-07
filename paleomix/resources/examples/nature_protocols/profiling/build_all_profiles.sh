#!/bin/bash
#
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
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
set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

{
## Pi1845A
echo Pi1845A.id.CGCTAT ../alignment/Pi1845A/reads/Pi1845A/Pi1845A_id_CGCTAT/*/reads.collapsed.bz2

## Pi1889
echo Pi1889.id.CTTGTA ../alignment/Pi1889/reads/Pi1889/Pi1889_id_CTTGTA/*/reads.collapsed.bz2
echo Pi1889.id.TAGCTT ../alignment/Pi1889/reads/Pi1889/Pi1889_id_TAGCTT/*/reads.collapsed.bz2
echo Pi1889.id.GGCTAC ../alignment/Pi1889/reads/Pi1889/Pi1889_id_GGCTAC/*/reads.collapsed.bz2

## M-0182896
echo M.0182896_UDG    ../alignment/M-0182896/reads/M-0182896/M-0182896_UDG/*/reads.collapsed.bz2
echo M.0182896_UDGa   ../alignment/M-0182896/reads/M-0182896/M-0182896_UDGa/*/reads.collapsed.bz2
echo M.0182896_UDGb   ../alignment/M-0182896/reads/M-0182896/M-0182896_UDGb/*/reads.collapsed.bz2
echo M.0182896_UDGc   ../alignment/M-0182896/reads/M-0182896/M-0182896_UDGc/*/reads.collapsed.bz2
echo M.0182896_NO_UDG ../alignment/M-0182896/reads/M-0182896/M-0182896_NO_UDG/*/reads.collapsed.bz2
} | xargs -P 16 -L1 ./build_profile.sh
