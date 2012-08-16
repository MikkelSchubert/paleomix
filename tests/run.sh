#!/bin/bash

# Ensure that timestamps are as expected
touch --date "2001-9-11 8:46" tests/data/timestamp_older
touch --date "2005-7-7 8:50"  tests/data/timestamp_younger

MODULES=$(ls *.py */*.py */*/*.py | grep -v "__init__" | sed -e 's#^#pypeline.#')

nosetests -d tests/ --with-coverage \
	--cover-erase \
	--cover-tests \
	--cover-inclusive \
	--cover-package=$(echo $MODULES | sed -e 's#\.py##g' -e's#/#.#g' -e's# #,#g')
