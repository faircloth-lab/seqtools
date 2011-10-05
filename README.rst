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

seqtools_ requires numpy_ (> 1.3).  After installing numpy, to install
seqtools_:

- from source::

    tar -xzvf ~/your/download/location/seqtools-*.tar.gz
    python setup.py install

- using easy_install::

    easy_install seqtools

- using pip::

    pip install seqtools


Tests
-----

While several of the functions and classes within seqtools_ are not (yet) 
fully tested, **ALL** of the *sequence-handling* classes, methods, and functions are
covered by unittests.

Running the tests as below, requires python-nose_.  After installing python-nose_, numpy_, 
and seqtools_ run the tests using:

>>> from seqtools import sequence
>>> sequence.test()


Alternatively, you can run:::

    python setup.py test

Or, after installing numpy_ and seqtools_ you can run the tests using:::

    python seqtools/sequence/tests/run.py


.. _Python: http://www.python.org/
.. _galaxy-tools: http://bitbucket.org/galaxy/galaxy-dist/src/
.. _numpy: http://numpy.scipy.org/
.. _lastz: http://www.bx.psu.edu/~rsharris/lastz/
.. _python-nose: http://code.google.com/p/python-nose/
.. _seqtools: https://github.com/faircloth-lab/seqtools/
