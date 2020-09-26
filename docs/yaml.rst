.. highlight:: YAML
.. _yaml_intro:

YAML usage in PALEOMIX
======================

`YAML`_ is a simple markup language adopted for use in configuration files by pipelines included in PALEOMIX. YAML was chosen because it is a plain-text format that is easy to read and write by hand. Since YAML files are plain-text, they may be edited using any standard text editors, with the following caveats:

* YAML exclusively uses spaces for indentation, not tabs; attempting to use tabs in YAML files will cause failures when the file is read by the pipelines.
* YAML is case-sensitive; an option such as `QualityOffset` is not the same as  `qualityoffset`.
* It is strongly recommended that all files be named using the `.yaml` file-extension; setting the extension helps ensure proper handling by editors that natively support the YAML format.

Only a subset of YAML features are actually used by PALEOMIX, which are described below. These include **mappings**, by which values are identified by names; **lists** of values; and **numbers**, **text-strings**, and **true** / **false** values, typically representing program options, file-paths, and the like. In addition, comments prefixed by the hash-sign (`#`) are frequently used to provide documentation.



Comments
--------

Comments are specified by prefixing unquoted text with the hash-sign (`#`); all comments are ignored, and have no effect on the operation of the program. Comments are used solely to document the YAML files used by the pipelines::

    # This is a comment; the next line contains both a value and a comment:
    123  # Comments may be placed on the same line as values.

For the purpose of the PALEOMIX reading this YAML code, the above is equivalent to the following YAML code::

    123

As noted above, this only applies to unquoted text, and the following is therefore not a comment, but rather a text-string::

    "# this is not a comment"

Comments are used in the following sections to provide context.


Numbers (integers and floats)
-----------------------------

Numbers in YAML file include whole numbers (integers) as well as real numbers (floating point numbers). Numbers are mostly used for program options, such as a minimum read length option, and involve whole numbers, but a few options do involve real numbers. Numbers may be written as follows::

    # This is an integer:
    123

    # This is a float:
    123.5

    # This is a float written using scientific notation:
    1.235e2


Truth-values (booleans)
-----------------------

Truth values (*true* and *false*) are frequently used to enable or disable options in PALEOMIX configuration files. Several synonyms are available which helps improve readability. More specifically, all of the following values are interpreted as *true* by the pipelines::

    true
    yes
    on

And similarly, the following values are all interpreted as *false*::

    false
    no
    off

Template files included with the pipelines mostly use `yes` and `no`, but either of the above corresponding values may be used. Note however that none of these values are quoted: If single or double-quotations were used, then these vales would be read as text rather than truth-values, as described next.


Text (strings)
--------------

Text, or strings, is the most commonly used type of value used in the PALEOMIX YAML files, as these are used to present both labels and values for options, including paths to files to use in an analysis::

    "Example"

    "This is a longer string"

    'This is also a string'

    "/path/to/my/files/reads.fastq"


For most part it is not necessary to use quotation marks, and the above could instead be written as follows::

    Example

    This is a longer string

    This is also a string

    /path/to/my/files/reads.fastq

However, it is important to make sure that values that are intended to be used strings are not misinterpreted as a different type of value. For example, without the quotation marks the following values would be interpreted as numbers or truth-values::

    "true"

    "20090212"

    "17e13"


Mappings
--------

Mappings associate a value with a label (key), and are used for the majority of options. A mapping is simply a label followed by a colon, and then the value associated with that label::

    MinimumQuality: 17

    EnableFoo: no

    NameOfTest: "test 17"

In PALEOMIX configuration files, labels are always strings, and are normally not quoted. However, in some cases, such as when using numerical labels in some contexts, it may be useful to quote the values:

    "A Label": on

    "12032016": "CPT"


Sections (mappings in mappings)
-------------------------------

In addition to mapping to a single value, a mapping may also itself contain one or more mappings::

    Top level:
      Second level: 'a value'
      Another value: true

Mappings can be nested any number of times, which is used in this manner to create sections and sub-sections in configuration files, grouping related options together::

    Options:
      Options for program:
        Option1: yes
        Option2: 17

      Another program:
        Option1: /path/to/file.fastq
        Option2: no

Note that the two mappings belonging to the `Option` mapping are both indented the same number of spaces, which is what allows the program to figure out which values belong to what label. It is therefore important to keep indentation consistent.

Lists of values
---------------

In some cases, it is possible to specify zero or more values with labels. This is accomplished using lists, which consist of values prefixed with a dash::

    Section:
      - First value
      - Second value
      - Third value

Note that the indentation of each item must be the same, similar to how indentation of sub-sections must be the same (see above).


Full example
------------

The following showcases basic structure of a YAML document, as used by the pipelines::

    # This is a comment; this line is completely ignored
    This is a section:
      This is a subsection:
        # This subsection contains 3 label / value pairs:
        First label: "First value"
        Second label: 2
        Third label: 3.14

      This is just another label: "Value!"

    This is a section containing a list:
      - The first item
      - The second item



.. _YAML: http://www.yaml.org
