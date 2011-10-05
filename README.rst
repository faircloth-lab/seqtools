Introduction
============

*seqtools* is a set of Python_ libraries for manipulating sequence data and 
generating sequence alignments.  Many of the functions and classes
within are derived from the excellent galaxy-tools_ libraries 
(`Blankenberg et al. <http://dx.doi.org/10.1093/bioinformatics/btq281>`_),
but modified to use numpy_ for storing quality values, etcetera.  There
are also functions for working with lastz_ and some general functions
for parsing arguments and validating input.

Generally speaking, this is a utility library that we use in other
programs on which we work.  Thus, it's commonly a dependency of many of
these packages.

Installation
------------

- from source
::
    tar -xzvf ~/your/download/location/seqtools.*.tar.gz

- using easy_install
::
    easy_install seqtools

- using pip
::
    pip install seqtools

Tests
-----

While several of the functions and classes within are not fully tested,
*ALL* of the sequence-handling classes, methods, and functions are
covered by unittests.  To run these tests after installation, do:

>>> from seqtools import sequence
>>> sequence.tests.run()

.. _Python: http://www.python.org/
.. _galaxy-tools: http://bitbucket.org/galaxy/galaxy-dist/src/
.. _numpy: http://numpy.scipy.org/
.. _lastz: http://www.bx.psu.edu/~rsharris/lastz/
