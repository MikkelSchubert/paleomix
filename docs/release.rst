Release checklist
-----------------

* Update changelog
* Update version in `paleomix/__init__.py` and in `docs/conf.py`
* Update version in paleomix_environment.yaml


Publish to PyPi
---------------

* git clone https://github.com/MikkelSchubert/paleomix.git
* cd paleomix
* tox
* git clean -fdx
* check-manifest
* python3 setup.py sdist
* twine upload -r testpypi dist/*
* twine upload dist/*


Publish to github
-----------------

* Tag release
* git push --tags
