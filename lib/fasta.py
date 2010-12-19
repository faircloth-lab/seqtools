#!/usr/bin/env python
# encoding: utf-8

"""
fasta.py

These files (fasta.py, fastq.py, transform.py, sequence.py) are originally part of the galaxy-dist 
package (http://bitbucket.org/galaxy/galaxy-dist/src/), and these particular classes, methods, and 
functions were created by Dan Blankenberg.  See Blankenberg et al. doi:  10.1093/bioinformatics/btq281.

I have modified the original source, adding methods for trimming and masking sequence reads, while
altering (slightly) how fasta and fastq files are read.  I've also added a class for reading
fasta+qual files.

Per the original license, the code was Copyright (c) 2005 Pennsylvania State University.

This file incorporates code covered by the following terms:

        Copyright (c) 2005 Pennsylvania State University

        Permission is hereby granted, free of charge, to any person obtaining
        a copy of this software and associated documentation files (the
        "Software"), to deal in the Software without restriction, including
        without limitation the rights to use, copy, modify, merge, publish,
        distribute, sublicense, and/or sell copies of the Software, and to
        permit persons to whom the Software is furnished to do so, subject to
        the following conditions:

        The above copyright notice and this permission notice shall be
        included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
        EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
        IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
        CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
        SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

import pdb
import copy
import numpy
from sequence import SequencingRead

class fastaSequence(SequencingRead):
    def __init__(self):
        SequencingRead.__init__(self)

class fastaReader(object):
    def __init__(self, fh):
        self.file = open(fh)
    
    def close( self ):
        """close file"""
        self.file.close()
    
    def next(self):
        """method for iterator"""
        line = self.file.readline()
        #remove header comment lines
        while line and line.startswith( '#' ): # pragma: no cover
            line = self.file.readline()
        if not line:
            raise StopIteration 
        assert line.startswith( '>' ), "FASTA headers must start with >"
        rval = fastaSequence()
        rval.identifier = line.strip()
        offset = self.file.tell()
        while True:
            line = self.file.readline()
            if not line or line.startswith( '>' ):
                if line:
                    self.file.seek( offset ) #this causes sequence id lines to be read twice, once to determine previous sequence end and again when getting actual sequence; can we cache this to prevent it from being re-read?
                return rval
            #454 qual test data that was used has decimal scores that don't have trailing spaces 
            #so we'll need to parse and build these sequences not based upon de facto standards
            #i.e. in a less than ideal fashion
            line = line.rstrip()
            if ' ' in rval.sequence or ' ' in line: # pragma: no cover
                rval.sequence = "%s%s " % ( rval.sequence, line )
            else:
                rval.sequence += line
            offset = self.file.tell()
    
    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()

class fastaQualityReader( object ):
    def __init__( self, fh, qual ):
        self.sequence_file = open(fh)
        self.quality_file = open(qual)
        
    def close( self ):
        """close file"""
        self.sequence_file.close()
        self.quality_file.close()
        
    def _seq( self, f, qual = False):
        """private method that is called to read both sequence and quality form f (input file)"""
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
    
    def next( self ):
        """method for iterator"""
        # read from sequence file until we hit a '>'
        sid, sequence = self._seq(self.sequence_file)
        qid, quality = self._seq(self.quality_file, True)
        assert sid == qid, "Seqeunce and quality identifiers are different"
        rval = fastaSequence()
        rval.identifier = sid
        rval.sequence = sequence
        rval.set_quality(quality)
        assert len(rval.sequence) == len(rval.quality), "Sequence and quality lengths are different"
        return rval

    def __iter__(self):
        """iterator"""
        while True:
            yield self.next()

class fastaWriter( object ):
    def __init__( self, fh, qual = False):
        self.sequence_file = open(fh, 'w')
        if qual:
            self.qual_file = open(qual, 'w')
        else:
            self.qual_file = qual
            
    def write( self, fasta, qual=False):
        """write fastaSequence objects to a file"""
        self.sequence_file.write(">{0}\n{1}\n".format(fasta.identifier[1:], fasta.sequence))
        if self.qual_file:
            self.qual_file.write(">{0}\n{1}\n".format(fasta.identifier[1:], fasta.get_quality_string()))
    
    def close( self ):
        """close output files"""
        self.sequence_file.close()
        if self.qual_file:
            self.qual_file.close()