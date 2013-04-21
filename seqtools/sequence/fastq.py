#!/usr/bin/env python
# encoding: utf-8

"""classes and methods for dealing with fastq-formatted files

Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

These files (fasta.py, fastq.py, transform.py, sequence.py) contain code and
ideas derived from galaxy-dist (http://bitbucket.org/galaxy/galaxy-dist/src/).
The original files were created by Dan Blankenberg.  See:

Blankenberg et al.  doi:  10.1093/bioinformatics/btq281

I have modified the original source, changing the way that quality scores are
stored and used (all quality scores are stored as numpy arrays in sanger-spec
integer (phred) values with methods to convert and display the array as a
standard quality string).  I've added additional methods to several of the
classes (e.g. trimming within sequence.SequenceRead - subclassed by both fasta
and fastq), and I altered methods and class methods to use numpy methods
within functions, for speed and where possible.  Finally, i slightly changed
how fasta and fastq files are read.

This file incorporates code covered by the following terms:

Copyright (c) 2005 Pennsylvania State University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import io
import math
import copy
import gzip
import numpy
import string
import transform
from os.path import splitext
from sequence import SequencingRead

import pdb

class FastqSequencingRead(SequencingRead):
    """Represents fastq sequences, attributes, and methods.  Default `type` is using Sanger-bases qualities"""
    format = 'sanger'
    ascii_min = 33
    ascii_max = 126
    quality_min = 0
    quality_max = 93
    score_system = 'phred'
    sequence_space = 'base'
    
    @classmethod
    def get_class_by_format(cls, format):
        """return the class of a fastq read and its associated attributes"""
        assert format in FASTQ_FORMATS.keys(), 'Unknown format type specified: %s' % format
        return FASTQ_FORMATS[ format ]
    
    @classmethod
    def convert_score_phred_to_solexa(cls, quality_array):
        """return an array of quality scores converted from the phred to solexa scoring system"""
        quality_array[numpy.where(quality_array <= 0)[0]] = -5
        return numpy.around(10.0 * numpy.log10(pow(10.0, quality_array/10.0) - 1.0)).astype('uint8')
    
    @classmethod
    def convert_score_solexa_to_phred(cls, quality_array):
        """return an array of quality scores converted from the solexa to phred scoring system"""
        return numpy.around(10.0 * numpy.log10(pow(10.0, quality_array/10.0) + 1.0)).astype('uint8')
    
    @classmethod
    def restrict_scores_to_valid_range(cls, quality_array):
        """return an array of scored limited to the range of the (format) class maximum"""
        return quality_array.clip(max=cls.quality_max)
    
    def set_quality(self, quality_string):
        """return ASCII and QUAL strings parsed as numpy.array([], dtype=uint8)"""
        quality_string = quality_string.strip()
        if ' ' in quality_string or '\t' in quality_string:
            self.quality = self._get_qual_array_from_decimal(quality_string)
        elif quality_string:
            temp_array = numpy.fromstring(quality_string, dtype = 'uint8')
            assert self.is_valid_format(temp_array), "Ascii quality string is not valid, given format."
            self.quality = self._get_qual_array_from_ascii(temp_array)
            self.assert_sequence_quality_lengths()
        else:
            self.quality = numpy.array([None])
    
    def _get_qual_array_from_ascii(self, quality_array):
        """PRIVATE: help parse an ASCII quality string and return an unsigned integer array of quality values"""
        return (quality_array - self.ascii_min + self.quality_min).astype('uint8')
    
    def get_quality_string(self, qtype = 'decimal'):
        """return the quality string of a given array in either decimal or ascii format"""
        if self.quality.any() and qtype == 'decimal':
            return ' '.join(self.quality.astype('|S2'))
        elif self.quality.any() and qtype == 'ascii':
            # numpy can return null bytes when there on non-letter characters in an array (e.g. 
            # "@".  strip those, return the string.)
            return (self.quality + self.ascii_min - (self.quality_min)).astype('uint8').tostring()
        else:
            raise ValueError, "Sequence has no quality values"

    def convert_read_to_format(self, format):
        """return an object converted from current format to requested format
        
        Convert an object from it's current format (self.format) to the
        requested format, which must be of the following type and scoring
        system: fasta, solexa, illumina.
        
        """
        assert format in FASTQ_FORMATS, 'Unknown format type specified: %s' % format
        new_class = FASTQ_FORMATS[ format ]
        new_read = new_class()
        new_read.identifier = self.identifier
        # we've removed color space
        new_read.sequence = self.sequence
        new_read.description = self.description
        if self.score_system == 'phred' and new_read.score_system == 'solexa':
            temp = self.convert_score_phred_to_solexa(self.quality)
        elif self.score_system == 'solexa' and new_read.score_system == 'phred':
            temp = self.convert_score_solexa_to_phred(self.quality)
        else:
            temp = self.quality
        new_read.quality = new_class.restrict_scores_to_valid_range(temp)
        return new_read
    
    def get_sequence(self):
        """return the string sequence of a read"""
        return self.sequence
    
    def is_valid_format(self, quality_string):
        """return BOOL indicating that ascii character values are within min-max for class/format"""
        if ((self.ascii_min <= quality_string) & (quality_string <= self.ascii_max)).all():
            return True
        else:
            return False
    
    def is_valid_quality_length( self ):
        """return BOOL indicating that sequence contains legitimate base pairs"""
        return len(self.quality) == len(self.sequence)
    
    def assert_sequence_quality_lengths(self):
        """ensure sequence length = quality length"""
        q,s = len(self.quality), len(self.sequence)
        assert s == q, "Invalid FASTQ file: quality score length ({0}) does not match sequence length ({1})".format(q,s)


class FastqSangerRead(FastqSequencingRead):
    """Represents valid attributed for sanger format sequence reads"""
    format = 'sanger'
    ascii_min = 33
    ascii_max = 126
    quality_min = 0
    quality_max = 93
    score_system = 'phred'
    sequence_space = 'base'

class FastqIlluminaRead(FastqSequencingRead):
    """Represents valid attributed for illumina format sequence reads"""
    format = 'illumina'
    ascii_min = 64
    ascii_max = 126
    quality_min = 0
    quality_max = 62
    score_system = 'phred'
    sequence_space = 'base'

class FastqSolexaRead(FastqSequencingRead):
    """Represents valid attributed for solexa format sequence reads"""
    format = 'solexa'
    ascii_min = 59
    ascii_max = 126
    quality_min = -5
    quality_max = 62
    score_system = 'solexa'
    sequence_space = 'base'

FASTQ_FORMATS = {}
for format in [ FastqIlluminaRead, FastqSolexaRead, FastqSangerRead ]:
    FASTQ_FORMATS[ format.format ] = format


class FastqSummarizer():
    """Represents a class for deriving summary data across a number of reads.
    Renamed from FastaAggregator in the original code
    
    >>> fastqs = fastq.FastqReader('../tests/test-data/sequence_for_summarizer.fastq')
    >>> fastq_stats = fastq.FastqSummarizer()
    >>> for read in fastqs:
    ...    fastq_stats.add(read)
    >>> # get summary values for a column - column[0] is the first base of all reads, 
    >>> # so we are getting a summary of the first base positions across all reads in 
    >>> # fastqs
    >>> summary_stats = fastq_stats.get_summary_for_column(0)
    
    
    """
    def __init__(self):
        self.lengths = {} #counts of seqs by read len
        self.nucleotides = ('A','C','G','T','N')
        self.positions = (0,1,2,3,4)
        self.nucleotide_positions = dict(zip(self.nucleotides, self.positions))
        self.bases = numpy.array([None]) # will be array to hold base counts
        self.qualities = numpy.array([None]) # we be array to hold quality stack
    
    def _count_bases(self, sequence):
        """PRIVATE: Return an array of base counts by column across all sequences
        aggregated.
        
        The array returned looks something like the following, when several 
        (n=2) sequences are aggregated:
        
            array([[ 0.,  1.,  0.,  1.,  0.],
                   [ 0.,  0.,  0.,  2.,  0.],
                   [ 0.,  0.,  0.,  2.,  0.],
                   ...
                   ])
        
        where each row represents a base position - thus row[0] is the count
        for sequence position 0.  The values within the vector give the counts
        of each possible base where the index is given in self.nucleotide_positions
        (in __init__):
        
            self.nucleotide_positions = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
        
        Using numpy arrays lets us quickly sum across columns or rows, making
        summary stats easier.
        
        """
        for k, base in enumerate(sequence):
            if numpy.ndim(self.bases) == 1 and k == 0:
                self.bases[self.nucleotide_positions[base]] += 1
            elif (numpy.ndim(self.bases) == 1) or (numpy.ndim(self.bases) > 1 and k == len(self.bases)):
                self.bases = numpy.vstack((self.bases, numpy.zeros(5)))
                self.bases[k][self.nucleotide_positions[base]] += 1
            else:
                self.bases[k][self.nucleotide_positions[base]] += 1
    
    def _nans(self, shape, dtype=numpy.float32):
        """PRIVATE helper method to return an array filled with NaNs"""
        a = numpy.empty(shape, dtype)
        a.fill(numpy.nan)
        return a
    
    def _array_qualities(self, quality):
        """PRIVATE: Return stacked quality arrays for reads into an array, buffering 
        read lengths to maximum encountered by appending numpy.NaNs to shorter
        reads.
        
        The returned array looks something like the following when several
        (n=3) sequences are aggregated:
        
            array([[ 17.,  17.,  17.,  17.,  17.,  17.,  17.],
                   [ 17.,  17.,  17.,  17.,  17.,  17.,  17.],
                   [ 17.,  17.,  17.,  17.,  17.,  17.,  17..]], dtype=float32)
        
        where each row of values represents the quality array of a sequence.
        Thus, we can min/max "down" columns to get the min/max over a particular
        base, across reads to get the per-read min/max, or over the entire
        array to get the overall min/max.
        
        """
        if not quality.any():
            raise ValueError, "There are no quality values"
        # deal with freaking switcheroo in numpy.array([]).shape
        if self.qualities.ndim == 1:
            diff = len(quality) - self.qualities.shape[0]
        else:
            diff = len(quality) - self.qualities.shape[1]
        if diff > 0:
            # hstack the necessary array of NaNs (e.g. diff x self.qualities 
            # length) - this is a little of a pain, but ensures we don't add
            # zeros to the array when the value is actually Null/None.
            self.qualities = numpy.hstack((self.qualities, 
                        self._nans([self.qualities.shape[0], diff])))
        # vstack quality
        self.qualities = numpy.vstack((self.qualities, quality))
    
    def add(self, read):
        """Add a read to the fastq aggregator
        
        In essence, we are building a "stack" of read data across all reads we 
        `add`. In order to generate summary statistics across some number of 
        reads, N, we need to add those reads to the aggregator, and then we 
        can compute our summary values across the reads that we add.
        
        >>> fastqs = fastq.FastqReader('../tests/test-data/sequence_for_summarizer.fastq')
        >>> fastq_stats = fastq.FastqSummarizer()
        >>> for read in fastqs:
        ...    fastq_stats.add(read)
        
        """
        # do a quick basecount
        if not self.bases.any():
            self.bases = numpy.zeros(5, dtype = numpy.float32)
        self._count_bases(read.sequence)
        # stack the quality arrays
        if not self.qualities.any():
            self.qualities = read.quality
        else:
            self._array_qualities(read.quality)
        # get dict of read lengths
        l = len(read)
        self.lengths[l] = self.lengths.get(l,0) + 1
    
    def get_decimal_range(self):
        """Return a tuple giving the (min, max) quality values for aggregated
        reads"""
        return (numpy.nanmin(self.qualities), numpy.nanmax(self.qualities))
    
    def get_length_counts(self):
        """Return a dictionary of read length counts"""
        return self.lengths
    
    def get_max_read_length(self):
        """Return the maximum read length"""
        return max(self.lengths.keys())
    
    def get_min_read_length(self):
        """Return the minumum read length"""
        return min(self.lengths.keys())
    
    def get_read_count_for_column(self, column):
        """Return the count of reads with a base in column"""
        if self.qualities.ndim > 1 and column >= self.qualities.shape[1]:
            return 0
        elif self.qualities.ndim == 1 and column >= self.qualities.shape[0]:
            return 0
        else:
            return sum(numpy.isfinite(self.qualities[:,column]))
    
    def get_read_count(self):
        """Return the total count of reads"""
        return self.get_read_count_for_column(0)
    
    def get_base_counts_for_column(self, column):
        """Return a dictionary of counts for each base in a column"""
        return dict(zip(self.nucleotides, self.bases[column]))
    
    def get_all_scores_for_column(self, column):
        """Return an array of quality scores for a column, including any NaNs"""
        return self.qualities[:, column]
    
    def get_finite_scores_for_column(self, column):
        """Return an array of quality scores for a column, excluding any NaNs"""
        return self.qualities[:, column][numpy.isfinite(self.qualities[:,column])]
    
    def get_nan_count_for_column(self, column):
        """Return the count of NaNs for a column"""
        return sum(numpy.isnan(self.qualities[:, column]))
    
    def get_quality_min_for_column(self, column):
        """Return the minimum quality score for a column (always exludes NaNs)"""
        return numpy.nanmin(self.qualities[:, column])
    
    def get_quality_max_for_column(self, column):
        """Return the maximum quality score for a column (always exludes NaNs)"""
        return numpy.nanmax(self.qualities[:, column])
    
    def get_quality_average_for_column(self, column):
        """Return the average quality score for a column (always excludes NaNs)"""
        return numpy.mean(self.get_finite_scores_for_column(column), dtype=numpy.float64)
    
    def get_quality_std_deviation_for_column(self, column):
        """Return the standard deviation of the quality for a column (always excludes NaNs)"""
        return numpy.std(self.get_finite_scores_for_column(column), ddof = 1, dtype=numpy.float64)
    
    def get_summary_for_column(self, column):
        reads = self.get_read_count_for_column(column)
        q_nans = self.get_nan_count_for_column(column)
        q_mn = self.get_quality_min_for_column(column)
        q_mx = self.get_quality_max_for_column(column)
        q_avg = self.get_quality_average_for_column(column)
        q_deviation = self.get_quality_std_deviation_for_column(column)
        return {'reads':reads, 'nans':q_nans, 'quality_min':q_mn, 
            'quality_max':q_mx, 'quality_avg':q_avg, 
            'quality_stddev':q_deviation, 'column':column
            }

class FastqReader():
    """Represents an iterator over fastaq sequences from a file"""
    def __init__(self, fastq_file, format = 'sanger'):
        #pdb.set_trace()
        if splitext(fastq_file)[1] == '.gz':
            self.file = io.BufferedReader(gzip.open(fastq_file, 'rb'))
        else:
            self.file = open(fastq_file)
        self.format = format
        
    def close(self):
        """close files"""
        return self.file.close()
        
    def next(self):
        """read next fastq sequence"""
        while True:
            fastq_header = self.file.readline()
            if not fastq_header:
                raise StopIteration
            fastq_header = fastq_header.rstrip('\n\r')
            #remove empty lines, apparently extra new lines at end of file is common?
            if fastq_header:
                break
        assert fastq_header.startswith('@'), 'Invalid fastq header: %s' % fastq_header
        rval = FastqSequencingRead.get_class_by_format(self.format)()        
        rval.identifier = fastq_header
        quality_string = ''
        while True:
            line = self.file.readline()
            if not line:
                raise Exception('Invalid FASTQ file: could not parse second instance of sequence identifier.')
            line = line.rstrip('\n\r')
            if line.startswith('+') and (len(line) == 1 or line[1:].startswith(fastq_header[1:])):
                rval.description = line
                break
            rval.append_sequence(line)
        while len(quality_string) < len(rval.sequence):
            line = self.file.readline()
            if not line:
                break
            quality_string = ''.join([quality_string, line.strip()])
        rval.set_quality(quality_string)
        rval.assert_sequence_quality_lengths()
        return rval
    
    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()

class FastqWriter():
    """Write fastq objects to a file"""
    def __init__(self, fastq_file, format = None):
        """set the sequence output file"""
        if splitext(fastq_file)[1] == '.gz':
            self.sequence_file = gzip.open(fastq_file, 'wb')
        else:
            self.sequence_file = open(fastq_file, 'w')
        self.format = format
    
    def write(self, fastq, qtype='ascii'):
        """write fastaSequence objects to a file"""
        if self.format:
            fastq = fastq.convert_read_to_format(self.format)
        self.sequence_file.write("{0}\n{1}\n+\n{2}\n".format(fastq.identifier, fastq.sequence, fastq.get_quality_string(qtype=qtype)))
    
    def close( self ):
        """close output files"""
        return self.sequence_file.close()

class FasterFastqReader(FastqReader):
    """A faster fastq reader - when you don't need quality arrays.

    Compared to the regular FastqWriter, this version is about 400 %
    faster.  Exection times for FasterFastqReader versus the regular
    FastqReader are below.  Sample sequence set was 100,000 lines.

    5.71s user 0.10s system 99% cpu 5.812 total
    22.96s user 0.10s system 99% cpu 23.075 total

    """
    def __init__(self, fastq_file, format = None):
        FastqReader.__init__(self, fastq_file, format = format)

    def next(self):
        """read next fastq sequence"""
        while True:
            fastq_header = self.file.readline()
            if not fastq_header:
                raise StopIteration
            fastq_header = fastq_header.rstrip('\n\r')
            #remove empty lines, apparently extra new lines at end of file is common?
            if fastq_header:
                break
        assert fastq_header.startswith('@'), 'Invalid fastq header: %s' % fastq_header
        quality_string = ''
        sequence = ''
        while True:
            line = self.file.readline()
            if not line:
                raise Exception('Invalid FASTQ file: could not parse second instance of sequence identifier.')
            line = line.rstrip('\n\r')
            if line.startswith('+') and (len(line) == 1 or line[1:].startswith(fastq_header[1:])):
                description = line
                break
            sequence = ''.join([sequence, line])
        while len(quality_string) < len(sequence):
            line = self.file.readline()
            if not line:
                break
            quality_string = ''.join([quality_string, line.strip()])
        return (fastq_header, description, sequence, quality_string)

class FasterFastqWriter(FastqWriter):
    """Write fastq objects to a file. NO conversion - just write
    the fastq out.  Take input form FasterFastqReader.
    """
    def __init__(self, fastq_file, format = None):
        FastqWriter.__init__(self, fastq_file, format = None)
    
    def write(self, fastq, qtype='ascii'):
        """write fastaSequence objects to a file"""
        #pdb.set_trace()
        identifier, description, sequence, quality = fastq
        self.sequence_file.write("{0}\n{1}\n{2}\n{3}\n".format(
                identifier,
                sequence,
                description,
                quality
            )
        )
