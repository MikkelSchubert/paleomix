Release checklist
-----------------

* Update changelog
* Update version in `paleomix/__init__.py` and in `docs/conf.py`
* Update version numbers referenced in `installation.rst`
* Update version in paleomix_environment.yaml


Publish to PyPi
---------------

* git clone https://github.com/MikkelSchubert/paleomix.git
* cd paleomix
* tox
* git clean -fdx
* check-manifest
* python3 setup.py sdist
* twine check dist/*
* twine upload -r testpypi dist/*
* pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple paleomix
* twine upload dist/*


Publish to github
-----------------

* Tag release
* git push --tags


Misc
----

* Prune list of versions listed on https://paleomix.readthedocs.io/
