#!/usr/bin/env python
# encoding: utf-8

"""classes and methods for dealing with fasta-formatted files

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

import copy
import numpy
import sequence

#import pdb

class FastaSequence(sequence.SequencingRead):
    """Represents fasta and fasta + qual sequences, attributes, and methods
    
    Can be used directly, but generally called from FastaReader or 
    fastaQualReader:
    
    >>> import fasta
    >>> seq = fasta.FastaSequence()
    >>> seq.identifier = 'chr5_6255117_6255601_0:0:0_1:0:0_13'
    >>> seq.sequence = 'CTTGGATCAGATGAAAATGCAGC'
    >>> q = "5 5 17 17 3 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"
    >>> seq.set_quality(q)
    >>> seq.sequence
    'CTTGGATCAGATGAAAATGCAGC'
    >>> seq.quality
    array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
           17, 17, 17, 17, 17, 17], dtype=uint8)
    """
    def __init__(self):
        sequence.SequencingRead.__init__(self)

class FastaReader():
    """Represents an iterator over fasta sequences from a file
    
    >>> import fasta
    >>> seq = '../tests/test-data/sequence.fasta'
    >>> fasta = fasta.FastaReader(seq)
    >>> seq = fasta.next()
    >>> seq.sequence
    'AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG'
    >>> fasta.close()
    
    """
    def __init__(self, fasta_file):
        """set the fasta file attribute"""
        self.sequence_file = open(fasta_file)
    
    def close( self ):
        """close file"""
        self.sequence_file.close()
    
    def next(self):
        """read next fasta sequence"""
        line = self.sequence_file.readline()
        #remove header comment lines
        while line and line.startswith( '#' ): # pragma: no cover
            line = self.sequence_file.readline()
        if not line:
            raise StopIteration 
        assert line.startswith( '>' ), "FASTA headers must start with `>`"
        rval = FastaSequence()
        rval.identifier = line.strip()
        offset = self.sequence_file.tell()
        while True:
            line = self.sequence_file.readline()
            if not line or line.startswith( '>' ):
                if line:
                    self.sequence_file.seek(offset)
                return rval
            line = line.rstrip()
            if ' ' in rval.sequence or ' ' in line: # pragma: no cover
                rval.sequence = "%s%s " % (rval.sequence, line)
            else:
                rval.sequence += line
            offset = self.sequence_file.tell()
    
    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()

class FastaQualityReader():
    """Represents an iterator over fasta+qual sequences from a file
    
    >>> import fasta
    >>> seq = '../tests/test-data/sequence.fasta'
    >>> qual = '../tests/test-data/sequence.qual'
    >>> fastaqual = fasta.FastaQualityReader(seq, qual)
    >>> seq = fastaqual.next()
    >>> seq.sequence
    'AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG'
    >>> seq.quality
    array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
           40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
           40, 40, 40, 40], dtype=uint8)
    >>> fastaqual.close()
    
    """
    def __init__( self, fasta_file, quality_file):
        """set the fasta file and qual attribute"""
        self.sequence_file = open(fasta_file)
        self.quality_file = open(quality_file)
        
    def close( self ):
        """close files"""
        self.sequence_file.close()
        self.quality_file.close()
        
    def _seq( self, f, qual = False):
        """PRIVATE: read both sequence and quality from input file"""
        line = f.readline()
        sequence = ''
        #remove header comment lines
        while line and line.startswith( '#' ): # pragma: no cover
            line = f.readline()
        if not line:
            raise StopIteration 
        assert line.startswith( '>' ), "FASTA headers must start with >"
        identifier = line.strip()
        offset = f.tell()
        while True:
            line = f.readline()
            if not line or line.startswith( '>' ):
                if line:
                    f.seek( offset )
                break
            #454 qual test data that was used has decimal scores that don't have trailing spaces 
            #so we'll need to parse and build these sequences not based upon de facto standards
            #i.e. in a less than ideal fashion
            if not qual:
                line = line.rstrip()
                assert ' ' not in line, "FASTA sequence contains spaces"
            if qual:
                # make sure we deal with corner case where qual line stops with `40\n` and starts
                # with `40`, giving us `4040` instead of `40 40`
                line = line.replace('\n', ' ')
            sequence += line
            offset = f.tell()
        if qual:
            # deal with any corner cases where we have extra spaces in the qual line
            sequence_list = sequence.rstrip().split(' ')
            while '' in sequence_list:
                sequence_list.remove('') # pragma: no cover
            # cast the list of strings to a list of ints
            sequence = ' '.join(sequence_list)
        return identifier, sequence
    
    def next(self):
        """read next fasta, qual sequence"""
        # read from sequence file until we hit a '>'
        sid, sequence = self._seq(self.sequence_file)
        qid, quality = self._seq(self.quality_file, True)
        assert sid == qid, "Seqeunce and quality identifiers are different"
        rval = FastaSequence()
        rval.identifier = sid
        rval.sequence = sequence
        rval.set_quality(quality)
        assert len(rval.sequence) == len(rval.quality), "Sequence and quality lengths are different"
        return rval

    def __iter__(self):
        """iterator"""
        while True:
            yield self.next()

class FastaWriter():
    """Write fasta and fasta+qual objects to a file
    
    >>> import fasta
    >>> seq = '../tests/test-data/sequence.fasta'
    >>> qual = '../tests/test-data/sequence.qual'
    >>> seq = fasta.FastaQualityReader(seq, qual)
    >>> outf = fasta.FastaWriter('../tests/test-data/test_write.fasta')
    >>> for s in seq:
    ...     outf.write(s)
    >>> outf.close()
    
    """
    def __init__(self, fasta_file, quality_file = False):
        """set the sequence and quality output files"""
        self.sequence_file = open(fasta_file, 'w')
        if quality_file:
            self.qual_file = open(quality_file, 'w')
        else:
            self.qual_file = quality_file
            
    def write(self, fasta, qual=False):
        """write FastaSequence objects to a file"""
        self.sequence_file.write(">{0}\n{1}\n".format(fasta.identifier.lstrip('>'), fasta.sequence))
        if self.qual_file:
            self.qual_file.write(">{0}\n{1}\n".format(fasta.identifier.lstrip('>'), fasta.get_quality_string()))
    
    def close(self):
        """close output files"""
        self.sequence_file.close()
        if self.qual_file:
            self.qual_file.close()

if __name__ == "__main__":
    import doctest # pragma: no cover
    doctest.testmod() # pragma: no cover
