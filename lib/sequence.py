#!/usr/bin/env python
# encoding: utf-8

"""
sequence.py

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

import numpy
import string
import transform
from copy import deepcopy

class SequencingRead(object):
    valid_sequence_list = string.letters
    
    def __init__( self ):
        self.identifier = None
        self.sequence = '' #holds raw sequence string: no whitespace
        self.description = None
        self.quality = '' #holds raw quality string: no whitespace, unless this contains decimal scores
        self.qual_array = numpy.array([None])
        self.trimming = False
        
    def __len__(self): # pragma: no cover
        return len(self.sequence)
    
    def __str__(self): # pragma: no cover
        return "%s\n%s\n%s\n%s\n" % (self.identifier, self.sequence, self.description, self.quality)
    
    def set_quality(self, quality_string):
        """parse ASCII and QUAL strings into and unsigned integer array that becomes the quality 
        attribute"""
        assert (' ' in quality_string or '\t' in quality_string), "Decimal quality values always contain spaces."
        self.quality = self._get_qual_array_from_decimal(quality_string)

    def get_quality_string(self):
        """return array as string"""
        return ' '.join(self.quality.astype('|S2'))
    
    def _get_qual_array_from_decimal(self, quality_string):
        """PRIVATE parse quality values 10 20 20 30 to array"""
        return numpy.array(quality_string.strip().split()).astype('uint8')
    
    def append_sequence(self, sequence):
        """append sequence"""
        self.sequence = ''.join([self.sequence, sequence.rstrip()])
    
    def append_quality(self, quality_string):
        """append quality"""
        assert ' ' in quality_string or '\t' in quality_string, "Decimal quality values always contain spaces."
        temp = self._get_qual_array_from_decimal(quality_string)
        self.quality = numpy.concatenate((self.quality,temp), axis = 0)
    
    def is_DNA(self):
        """sequence is not RNA"""
        return 'u' not in self.sequence.lower()
    
    def is_RNA(self):
        """sequence is not DNA"""
        return 'u' in self.sequence.lower()
    
    def clone(self):
        """create unlinked copy of object"""
        return deepcopy(self)
        
    def snapshot(self):
        """snapshot sequence and quality"""
        self.sequence_snapshot = deepcopy(self.sequence)
        self.quality_snapshot = deepcopy(self.quality)
    
    def reverse(self, clone = True):
        """reverse sequence strand"""
        if clone:
            rval = self.clone()
        else:
            rval = self
        rval.sequence = transform.reverse(self.sequence)
        rval.quality = rval.quality[::-1]
        return rval
    
    def complement(self, clone = True):
        """complement sequence strand"""
        if clone:
            rval = self.clone()
        else:
            rval = self
        if rval.is_DNA():
            rval.sequence = transform.DNA_complement(rval.sequence)
        else:
            rval.sequence = transform.RNA_complement(rval.sequence)
        return rval
    
    def reverse_complement(self, clone = True):
        """reverse complement sequence strand"""
        #need to reverse first, then complement
        rval = self.reverse(clone = clone)
        return rval.complement(clone = False) #already working with a clone if requested
    
    def sequence_as_DNA(self, clone = True):
        """convert RNA to DNA"""
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.sequence = transform.to_DNA(rval.sequence)
        return rval
        
    def sequence_as_RNA(self, clone = True):
        """convert DNA to RNA"""
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.sequence = transform.to_RNA(rval.sequence)
        return rval
    
    def _check(self):
        """PRIVATE:  determine which quality scores fall below self.min_qual"""
        comparison = self.quality < self.min_qual
        if True not in comparison:
            return False, None
        else:
            return True, comparison

    def mask_and_trim(self, min_qual):
        """mask min_qual bases to N.  trim min_qual bases from ends"""
        assert self.quality.any(), "There is no quality information to use for trimming"
        self.min_qual = min_qual
        trim, comparison = self._check()
        if trim:
            self.snapshot()
            sequence_array = numpy.array(list(self.sequence))
            # get values where the qual is low and replace them w/ N
            sequence_array[numpy.where(comparison)[0]] = 'N'
            # trim and adjust self.sequence
            self.sequence = sequence_array.tostring().lstrip('N').rstrip('N')
            # trim and adjust self.quality
            temp = deepcopy(self.quality)
            temp[numpy.where(comparison)[0]] = 0
            self.quality[numpy.where(comparison)[0]] = 0
            self.quality = numpy.trim_zeros(self.quality)
            # ensure the trimming didn't remove something inside the read
            assert len(self.sequence) == len(self.quality)
            self.trimming = 'mt'

    def mask(self, min_qual, ambiguous = 'N'):
        """mask min_qual bases to N.  no trimming"""
        assert self.quality.any(), "There is no quality information to use for trimming"
        self.min_qual = min_qual
        trim, comparison = self._check()
        if trim:
            self.snapshot()
            sequence_array = numpy.array(list(self.sequence))
            sequence_array[numpy.where(comparison)[0]] = ambiguous
            self.quality[numpy.where(comparison)[0]] = 0
            self.sequence = sequence_array.tostring()
            self.trimming = 'm'

    def trim(self, min_qual):
        """trim min_qual bases from ends.  no inner sequence masking"""
        assert self.quality.any(), "There is no quality information to use for trimming"
        self.min_qual = min_qual
        trim, comparison = self._check()
        if trim:
            self.snapshot()
            # use temp array so we don't change inner quality values (we are *not* masking)
            temp = deepcopy(self.quality)
            temp[numpy.where(comparison)[0]] = 0
            l = len(self.quality) - len(numpy.trim_zeros(temp,'f'))
            r = len(self.quality) - len(numpy.trim_zeros(temp,'b'))
            self.quality = self.quality[l:len(self.quality) - r]
            self.sequence = self.sequence[l:len(self.sequence) - r]
            self.trimming = 't'
    
    def slice(self, left_column_offset, right_column_offset, clone = True):
        """slice a sequence between left and right, returning a new object, if requested"""
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.sequence = rval.sequence[left_column_offset:right_column_offset]
        if rval.quality.any():
            rval.quality = rval.quality[left_column_offset:right_column_offset]
        return rval
        
