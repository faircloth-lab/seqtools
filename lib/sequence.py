#!/usr/bin/env python
# encoding: utf-8

"""establishes the SequencingRead class and methods 

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

import numpy
import string
import transform
from copy import deepcopy

class SequencingRead():
    """A generic class encapuslating DNA sequence related methods and attributes"""
    valid_sequence_list = string.letters
    def __init__(self):
        """initialize the SequencingRead class"""
        self.identifier = None
        self.sequence = '' #holds raw sequence string: no whitespace
        self.description = None
        self.quality = numpy.array([None])
        self.trimming = False
        
    def __len__(self): # pragma: no cover
        """return length of self.sequence"""
        return len(self.sequence)
    
    def __str__(self): # pragma: no cover
        """return representation of the object"""
        return "%s\n%s\n%s\n%s\n" % (self.identifier, self.sequence, self.description, self.quality)
    
    def set_quality(self, quality_string):
        """set self.quality attribute
        
        parse ASCII and QUAL quality_string into an unsigned integer array
        that becomes the quality attribute
        
        """
        assert (' ' in quality_string or '\t' in quality_string), "Decimal quality values always contain spaces."
        self.quality = self._get_qual_array_from_decimal(quality_string)

    def get_quality_string(self):
        """return self.quality array as decimal string"""
        return ' '.join(self.quality.astype('|S2'))
    
    def _get_qual_array_from_decimal(self, quality_string):
        """PRIVATE: parse quality values 10 20 20 30 to array"""
        return numpy.array(quality_string.strip().split()).astype('uint8')
    
    def append_sequence(self, sequence):
        """append sequence to self.sequence attribute"""
        self.sequence = ''.join([self.sequence, sequence.rstrip()])
    
    def append_quality(self, quality_string):
        """append quality to self.quality attribute"""
        assert ' ' in quality_string or '\t' in quality_string, "Decimal quality values always contain spaces."
        temp = self._get_qual_array_from_decimal(quality_string)
        self.quality = numpy.concatenate((self.quality,temp), axis = 0)
    
    def is_DNA(self):
        """return BOOL indicating if sequence is not RNA"""
        return 'u' not in self.sequence.lower()
    
    def is_RNA(self):
        """return BOOL indicating if sequence is not DNA"""
        return 'u' in self.sequence.lower()
    
    def is_valid_sequence(self):
        """return BOOL indicating that sequence contains legitimate base pairs"""
        valid_bases = set(['a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'r', 'u', 't', 'v', 'y'])
        temp = set([i for i in self.sequence.lower()])
        return temp.issubset(valid_bases)
    
    def clone(self):
        """return a copy.deepcopy (clone) of this object
        
        we use this when we are going to make any potentially destructive
        changes to an object (e.g. trimming, masking, slicing).  The default 
        case is to *always* return a clone of the object rather than modifying
        the attributes of the object in-place.  This is different than
        self.snapshot(), in that snapshot just archives the initial state
        of sequence and quality attributes before manipulation.
        
        """
        return deepcopy(self)
        
    def snapshot(self):
        """archive untouched copies of self.sequence and self.quality
        
        archive untouched copies of self.sequence and self.quality as 
        self.sequence_snapshot and self.quality_snapshot. typically called 
        before we run some 'destructive' manipulation on the object in 
        question (e.g. trimming).  Different from self.clone() in that
        clone returns a copy of the entire object.
        
        """
        self.sequence_snapshot = deepcopy(self.sequence)
        self.quality_snapshot = deepcopy(self.quality)
    
    def reverse(self, clone = True):
        """return a new object with a reverses self.sequence strand
        
        call using my_object.reverse(clone=False) if you wish to reverse and
        return the current object.
        
        """
        if clone:
            rval = self.clone()
        else:
            rval = self
        rval.sequence = transform.reverse(self.sequence)
        rval.quality = rval.quality[::-1]
        return rval
    
    def complement(self, clone = True):
        """return a new object with a complemented self.sequence strand
        
        call using my_object.complement(clone=False) if you wish to complement
        and return the current object.
        
        """
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
        """return a new object with a reverse-complemented self.sequence strand
        
        call using my_object.reverse_complement(clone=False) if you wish to 
        complement and return the current object.
        
        """
        #need to reverse first, then complement
        rval = self.reverse(clone = clone)
        return rval.complement(clone = False) #already working with a clone if requested
    
    def sequence_as_DNA(self, clone = True):
        """return a new object with an RNA self.sequence strand as DNA
        
        call using my_object.sequence_as_DNA(clone=False) if you wish to 
        convert and return the current object.
        
        """
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.sequence = transform.to_DNA(rval.sequence)
        return rval
        
    def sequence_as_RNA(self, clone = True):
        """return a new object with a DNA self.sequence strand as RNA
        
        call using my_object.sequence_as_RNA(clone=False) if you wish to 
        convert and return the current object.
        
        """
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

    def _qual_check(self):
        """PRIVATE:  return BOOL indicating whether we have quality values"""
        assert self.quality.any(), "There is no quality information to use for trimming"
    
    def mask_and_trim(self, min_qual, clone = True):
        """return a new object with DNA masked and trimmed according to some min_qual
            
        call using my_object.mask_and_trim(clone=False) if you wish to 
        mask, trim, and return the current object.
        
        """
        self._qual_check()
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.min_qual = min_qual
        trim, comparison = rval._check()
        if trim:
            rval.snapshot()
            sequence_array = numpy.array(list(rval.sequence))
            # get values where the qual is low and replace them w/ N
            sequence_array[numpy.where(comparison)[0]] = 'N'
            # trim and adjust rval.sequence
            rval.sequence = sequence_array.tostring().lstrip('N').rstrip('N')
            # trim and adjust rval.quality
            temp = deepcopy(rval.quality)
            temp[numpy.where(comparison)[0]] = 0
            rval.quality[numpy.where(comparison)[0]] = 0
            rval.quality = numpy.trim_zeros(rval.quality)
            # ensure the trimming didn't remove something inside the read
            assert len(rval.sequence) == len(rval.quality)
            rval.trimming = 'mt'
        return rval

    def mask(self, min_qual, ambiguous = 'N', clone = True):
        """return a new object with DNA masked according to some min_qual
            
        call using my_object.mask(clone=False) if you wish to 
        mask, trim, and return the current object.  You may also specify
        the masking character for bases < min_qual using 
        my_object.mask(ambiguous = '?')
        
        """
        self._qual_check()
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.min_qual = min_qual
        trim, comparison = rval._check()
        if trim:
            rval.snapshot()
            sequence_array = numpy.array(list(rval.sequence))
            sequence_array[numpy.where(comparison)[0]] = ambiguous
            rval.quality[numpy.where(comparison)[0]] = 0
            rval.sequence = sequence_array.tostring()
            rval.trimming = 'm'
        return rval

    def trim(self, min_qual, clone = True):
        """return a new object with DNA trimmed according to some min_qual
            
        call using my_object.trim(clone=False) if you wish to 
        mask, trim, and return the current object.
        
        """
        self._qual_check()
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.min_qual = min_qual
        trim, comparison = rval._check()
        if trim:
            rval.snapshot()
            # use temp array so we don't change inner quality values (we are *not* masking)
            temp = deepcopy(rval.quality)
            temp[numpy.where(comparison)[0]] = 0
            l = len(rval.quality) - len(numpy.trim_zeros(temp,'f'))
            r = len(rval.quality) - len(numpy.trim_zeros(temp,'b'))
            rval.quality = rval.quality[l:len(rval.quality) - r]
            rval.sequence = rval.sequence[l:len(rval.sequence) - r]
            rval.trimming = 't'
        return rval
    
    def slice(self, left_column_offset, right_column_offset, clone = True):
        """return a new object with DNA sliced from left_column_offset to right_column_offset
            
        call using my_object.slice(10, 20, clone=False) if you wish to 
        mask, trim, and return the current object.
        
        """
        if clone: # pragma: no cover
            rval = self.clone()
        else: # pragma: no cover
            rval = self
        rval.sequence = rval.sequence[left_column_offset:right_column_offset]
        if rval.quality.any():
            rval.quality = rval.quality[left_column_offset:right_column_offset]
        return rval
        
