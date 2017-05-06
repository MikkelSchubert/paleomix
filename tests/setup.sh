#!/bin/bash

set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

if [ ! -e "paleomix" ];
then
    cd ..
fi

find paleomix -name '*.py' -type f \
    | sed -e's#\.py##' -e's#/#.#g' -e's#^#import #' \
    | grep -v "__init__\|yaml\|resources" \
    > tests/all_modules.py