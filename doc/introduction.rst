.. _introduction:

***************
Introduction
***************

.. _next_generation_sequencing:

Next-generation sequencing
===========================

'Next-generation sequencing' or 'next-generation sequencers' (very) generally refer to techniques yielding massive (order of magnitude) increases in the acquisition of DNA sequence data and concomitant reductions in the price-per-base of sequence produced.  The methodologies used by these different platforms vary (sometimes vastly), and you should read something along the lines of [NextGenReview]_ for additional, technical information.

Platforms
----------

Throughout the following and with respect to "next-generation sequencers", we are generally speaking of the Roche `454 FLX Series <http://454.com/products-solutions/system-benefits.asp>`_, Illumina's `Genome Analyser <http://www.illumina.com/systems/genome_analyzer_iix.ilmn>`_, and ABI's `SOLiD system <http://www3.appliedbiosystems.com/AB_Home/applicationstechnologies/SOLiDSystemSequencing/index.htm>`_.

There exists a third category of DNA sequencers now (2010) that we could refer to as next-next-generation sequencers, as these represent the third major iteration of sequencer design.  Included among this group are higher-throughput sequencers from `Pacific Biosciences <http://www.pacificbiosciences.com/>`_,  `Illumina  <http://www.illumina.com/systems/hiseq_2000.ilmn>`_ and lower-cost sequencers from `ion torrent <http://www.iontorrent.com>`_.

Long-read v. Short-read
------------------------

You may see next-generation sequencers referred to as long-read or short read platforms.  Generally speaking, the Roche 454 and Pacific Biosciences machines are considered "long-read" machine, because the DNA reads they return can vary from 300 to 1000+ bp in length, while the machines from Illumina and ABI are considered "short-read" machines - the reads returned by these instruments are typically from 50 - 100 bp in length. 

The Next-gen Sequencing 'Problem'
=================================

Next-generation DNA sequencers bring along vast increases in the amount of read data returned from a single run of a particular machine.  For instance, one run of the (first-gen) Applied Biosystems 3730xl capillary DNA sequencer would return an entire plate (96 wells) of sequence reads, each approximately 600 bases in length.  All together, that totals 56,700 (96 x 600) base pairs (57 kb) of sequence returned per run.

Next-generation sequencers, on the other hand, typically return anywhere from 400-600 million bases (0.5 gb or 400,000 kb or 7000x that of the 3730xl) for the Roche 454 FLX with Titanium chemistry to 150-200 billion bases (150-200 gb) per run on the recently announced Illumina HiSeq [`link <http://www.illumina.com/systems/hiseq_2000.ilmn>`_].

The count of reads returned is generally the attraction of these new machines - we can collect massive amounts of data from a single run.  However, many of these machines start with the premise that one is interested in collecting this count of reads from the same (or a handful) of individuals or organisms.  Thus, the data returned often provide high **depths** of sequencing coverage from but a few individuals.

Many biologists, while interested in depth of coverage, are not necessarily interested in the extreme depths of coverage necessary to type extremely rare variants of a particular gene/genome region of interest.  Rather, we would like to be able to take the **depth** of coverage provided by these machines spread it out across many individuals or a *population* of individuals - giving us **breadth** of coverage at the cost of **depth** of coverage.

.. _what_is_sequence_tagging:

What is sequence tagging
=========================

Sequence tagging [#f1]_ is a means of using synthetic pieces of DNA attached to DNA reads of interest to help identify those pieces of DNA after they have been sequenced.  In essence, they provide a reference, attached to each read, that allows us, later, to determine the individual/population from which the DNA read came.

.. [#f1] sequence tagging is also called:  a) molecular tagging; b) molecular barcoding; c) sequence barcoding; etc.