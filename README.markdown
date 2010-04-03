# linker-py

Copyright (c) 2009-2010 Brant C. Faircloth.  All rights reserved.

## Description

`linker-py` is a program for the parsing, error-correcting, and tracking of DNA reads identified using molecular sequence tags.  `linker-py` differs from previous software/methodological approaches because:

 * it uses fuzzy string matching along with sequence tags of Levenshtein distance >= X to identify and correct sequencing errors, which are somewhat common (1%) on certain next-generation sequencing platforms.  The Levenshtein distance differs from other implementations (i,e. Hamming distance) in that the distance represents the number of insertions, deletions and/or substitutions needed to get from one sequence of characters to another.
 
 * it can parse, error-correct, and track "traditional" sequence tagged reads in addition to hierarchically sequence-tagged reads, which vastly expand the number of sequence pools that may be mixed during any one next-generation sequencing run
 
 * it organizes sequence read data, by tag or other metadata, in a relational database for additional downstream processing
 
 * it can take advantage of multiple processors/cores for the parsing and error correcting of sequence reads to reduce overall processing time
 
 * it is flexible with respect to the linker sequences used

## Dependencies

* Python ( >2.6 or 2.5 with multiprocessing backport )
* MySQL ( tested with > version 5.x.x )
* MySQLdb (v. )
* Biopython ( > v1.54 )
* progress.py ( available in linker-py package or )
* a valid linker-py.conf file


## Installation

* copy the script to a location of your choice
* you may choose to `chmod 0755` the script after placing it in your $PATH

## Invocation

You run the program using:

    python --configuration=linker-py.conf`


## Extras

`linker-py` includes several additional 'helper' scripts to make your life easier:

 * `levenshtein.py` - this script will take a linker-py.conf file and evaluate the Levenshtein distance (or Hamming distance) between the tags on a per-section basis.
        
        python levenshtein.py --configuration=linker-py.conf --section=Linker

 * `getLinkerSequence..py` - this script will retrieve sequence from the database given a number of criteria
        
        python getLinkerSequence.py --configuration=linker-py.conf
