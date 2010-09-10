# linker-py

Copyright (c) 2009-2010 Brant C. Faircloth.
All rights reserved.

See `LICENSE.markdown` for standard, 2-clause BSD license.

## Description

`linker-py` is a program for the parsing, error-correcting, and tracking of DNA reads identified using molecular sequence tags.  `linker-py` differs from previous software/methodological approaches because:

 * it uses fuzzy string matching along with sequence tags of Levenshtein distance >= X to identify and correct sequencing errors, which are somewhat common (~1%; Roche, Inc. personal communication) on certain next-generation sequencing platforms.  The [Levenshtein][1] distance differs from other implementations (i,e. [Hamming][2] distance) in that the distance represents the number of insertions, deletions and/or substitutions needed to get from one sequence of characters to another.  For additional information see Suggested Readings, below.
 
 * it can parse, error-correct, and track "traditional" (e.g. [Meyer et al. 2007][3]) sequence tagged reads in addition to hierarchically sequence-tagged reads (e.g. [Glenn et al. 2010][8]).  Hierarchical tagging vastly expands the number of sequence pools that may be mixed during any single next-generation sequencing run
 
 * it organizes sequence read data, by tag or other metadata, in a relational database ([Mysql][9]) for additional downstream processing
 
 * it can take advantage of [multiprocessing][10] for the parsing and error correcting of sequence reads to reduce overall processing time
 
 * linkers and hierarchical combinations are flexible and specific in a easily edited [configuration file][11]

## Dependencies

* Python ( tested with 2.6.3 and 2.6.4; works with 2.5.x and multiprocessing backport )
* numpy ( tested with 1.3.0 )
* MySQL ( tested with > version 5.x.x )
* MySQLdb ( tested with 5.1.45 )
* Biopython ( tested with 1.53 )
* progress.py ( available in linker-py package )
* pylevenshtein ( for helper/levenshtein.py available from [http://code.google.com/p/pylevenshtein/](http://code.google.com/p/pylevenshtein/) )
* a valid linker-py.conf file ( see example )


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

 * `getLinkerSequence.py` - this script will retrieve sequence from the database given a number of criteria
        
        python getLinkerSequence.py --configuration=linker-py.conf
        
## Suggested Reading

 * Sequence tagging
    * Meyer M, Stenzel U, Myles S, Prüfer K, Hofreiter M (2007) Targeted 
    high-throughput sequencing of tagged nucleic acid samples.  Nucleic Acids 
    Research 35(15):e97.  [LINK][3]
    * Hamady M, Walker JJ, Harris JK, Gold NJ, Knight R (2008)
    Error-correcting barcoded primers for pyrosequencing hundreds of samples 
    in multiplex.  Nature Methods 5 (3):235-237. 
    [LINK][4]
    * Binladen J, Gilbert MT, Bollback JP, Panitz F, Bendixen C, Nielsen R, 
    Willerslev E (2009) The use of coded PCR primers enables high-throughput 
    sequencing of multiple homolog amplification products by 454 parallel 
    sequencing.  BMC Bioinformatics 10:362. 
    [LINK][5]
 * Levenshtein distance
     *  Levenshtein VI (1966). Binary codes capable of correcting deletions, 
     insertions, and reversals. Soviet Physics Doklady 10:707–10. [LINK][6]
 * Hamming distance
     * Hamming RW (1950) Error detecting and error correcting codes. Bell 
     System Technical Journal 26 (2):147–160. [LINK][7]
     
        
[1]:  http://en.wikipedia.org/wiki/Levenshtein_distance
[2]:  http://en.wikipedia.org/wiki/Hamming_distance
[3]:  http://dx.doi.org/10.1093/nar/gkm566
[4]:  http://dx.doi.org/10.1038/nmeth.1184
[5]:  http://dx.doi.org/10.1371/journal.pone.0000197
[6]:  http://sascha.geekheim.de/wp-content/uploads/2006/04/levenshtein.pdf
[7]:  http://www.caip.rutgers.edu/~bushnell/dsdwebsite/hamming.pdf
[8]:  http://www.uga.edu/
[9]:  http://www.mysql.com/
[10]: http://en.wikipedia.org/wiki/Multiprocessing
[11]: https://github.com/brantfaircloth/linker-py/blob/master/linker-py.conf