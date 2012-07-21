#!/bin/bash

MODULES=$(ls *.py */*.py */*/*.py | grep -v "__init__" | sed -e 's#^#pypeline.#')

nosetests -d tests/ --with-coverage \
	--cover-erase \
	--cover-tests \
	--cover-inclusive \
	--cover-package=$(echo $MODULES | sed -e 's#\.py##g' -e's#/#.#g' -e's# #,#g')
