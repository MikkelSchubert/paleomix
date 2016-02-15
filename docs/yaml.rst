.. highlight:: YAML
.. _yaml_intro:

YAML usage in PALEOMIX
======================

The format, `YAML`_, is a simple human-readable markup language in which the structure of the data is determined by its identation, and will look familiar to anyone who has experience with the `Python`_ programming language.

The following showcases basic structure of a YAML document, as used by the pipelines::

    # This is a comment; this line is completely ignored
    This is a section:
      This is a subsection:
        # This subsection contains 4 key / value pairs:
        First key: "First value"
        Second key: 2
        Third key: 3.4
        # The following key has no value assosiated!
        Fourth key:

    This is a section containing a list:
      - The first item
      - The second item



.. _Python: http://www.python.org/
.. _YAML: http://www.yaml.org
