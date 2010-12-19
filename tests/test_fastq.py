#!/usr/bin/env python
# encoding: utf-8

"""
test_fastq.py

Created by Brant Faircloth on 14 December 2010 22:19 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.


This file incorporates several snippets of code covered by the following terms:

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

import os
import numpy
import math
import unittest
import exceptions
from demuxipy.lib import fastq

import pdb

class TestFastqClassMethods(unittest.TestCase):
    def setUp(self):
        self.s = fastq.fastqSequencingRead()
        self.s.identifier = "@chr5_6255117_6255601_0:0:0_1:0:0_13/1"
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
        self.s.set_quality("&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$")
    
    def test_class_format(self):
        """[fastq] check for incorrect format"""
        self.assertRaises(AssertionError, self.s.get_class_by_format, 'higgs')
    
    def old_phred_to_solexa(self, score):
        """original classmethod in Dan Blankenberg's code"""
        if score <= 0: #can't take log10( 1 - 1 ); make <= 0 into -5
            return -5
        return int(round(10.0 * math.log10(math.pow(10.0, (float(score)/10.0)) - 1.0)))

    def old_convert_score_phred_to_solexa(self, decimal_score_list):
        """original classmethod in Dan Blankenberg's code"""
        return map(self.old_phred_to_solexa, decimal_score_list)

    def test_convert_score_phred_to_solexa(self):
        """[fastq] convert phred scoring to solexa scoring"""
        old = numpy.array(self.old_convert_score_phred_to_solexa(self.s.quality.tolist()))
        new = self.s.convert_score_phred_to_solexa(self.s.quality)
        assert (new == old).all()
    
    def old_solexa_to_phred(self, score):
        """original classmethod in Dan Blankenberg's code"""
        return int(round(10.0 * math.log10(math.pow(10.0, (float(score)/10.0)) + 1.0)))
    
    def old_convert_score_solexa_to_phred(self, decimal_score_list):
        """original classmethod in Dan Blankenberg's code"""
        return map(self.old_solexa_to_phred, decimal_score_list)
    
    def test_convert_score_solexa_to_phred(self):
        """[fastq] convert solexa scoring to phred scoring"""
        old = numpy.array(self.old_convert_score_solexa_to_phred(self.s.quality.tolist()))
        new = self.s.convert_score_solexa_to_phred(self.s.quality)
        assert (new == old).all()
    
    def old_restrict_score(self, score):
        """original classmethod in Dan Blankenberg's code"""
        return max(min(score, self.s.quality_max), self.s.quality_min)
    
    def old_restrict_scores_to_valid_range(self, decimal_score_list):
        """original classmethod in Dan Blankenberg's code"""
        return map(self.old_restrict_score, decimal_score_list)
    
    def test_restrict_scores_to_valid_range(self):
        """[fastq] ensure that scores do not exceed the maximum for the standard"""
        self.s.quality = numpy.array([100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                100, 100, 100, 90, 90, 80, 70])
        old = self.old_restrict_scores_to_valid_range(self.s.quality.tolist())
        new = self.s.restrict_scores_to_valid_range(self.s.quality)
        assert (new == old).all()
    
    def old_get_ascii_quality_scores(self, quality):
        """original ascii conversion in Dan Blankenberg's code"""
        quality = quality.rstrip() #decimal scores should have a trailing space
        if quality:
            try:
                return [ chr( int( val ) + self.s.ascii_min - self.s.quality_min ) for val in quality.split() ]
            except ValueError, e:
                raise ValueError( 'Error Parsing quality String. ASCII quality strings cannot contain spaces (%s): %s' % ( self.quality, e ) )
        else:
            return []
    
    def test_get_ascii_quality_scores(self):
        """[fastq] convert decimal scores to ascii scores"""
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
        q = '10 10 20 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 '
        self.s.set_quality(q)
        old = ''.join(self.old_get_ascii_quality_scores(q))
        new = self.s.get_quality_string('ascii')
        assert new == old
    
    def test_get_decimal_quality_scores(self):
        """[fastq] convert ascii scores to decimal scores"""
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
        self.s.set_quality('++5?IIIIIIIIIIIIIII')
        assert (self.s.quality == [10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]).all()
        self.s.set_quality('10 10 20 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 ')
        assert (self.s.quality == [10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]).all()
        self.s.set_quality(' ')
        assert self.s.quality == None
    
class TestFastqConversions(unittest.TestCase):
    
    def setUp(self):
        self.identifier = '@chr5_56045273_56045722_1:0:0_2:0:0_12/1'
        self.sequence = 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        self.sanger_quality = '&&22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$'
        self.illumina_quality = 'EEQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQEDC'
        self.solexa_quality = 'CCQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQCB@'
        self.expected_quality = numpy.array([ 5,  5, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])
        self.expected_solexa_quality = numpy.array([ 3,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                      17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                      17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                      17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                      17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                      17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  3,  2,  0])
    
    def test_for_wrong_format(self):
        """[fastq] disallow improper quality formats"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('illumina')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.illumina_quality)
        self.assertRaises(AssertionError, self.observed.convert_read_to_format, 'bob')

    def test_convert_ilumina_to_sanger(self):
        """[fastq] convert illumina to sanger"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('illumina')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.illumina_quality)
        # the array value of this should be equal to our "usual" array value
        assert (self.observed.quality == self.expected_quality).all()
        converted = self.observed.convert_read_to_format('sanger')
        assert converted.get_quality_string('ascii') == self.sanger_quality
    
    def test_convert_sanger_to_illumina(self):
        """[fastq] convert sanger to illumina"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('sanger')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.sanger_quality)
        assert (self.observed.quality == self.expected_quality).all()
        converted = self.observed.convert_read_to_format('illumina')
        assert converted.get_quality_string('ascii') == self.illumina_quality
    
    def test_convert_sanger_to_solexa(self):
        """[fastq] convert sanger to solexa"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('sanger')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.sanger_quality)
        converted = self.observed.convert_read_to_format('solexa')
        assert (converted.quality == self.expected_solexa_quality).all()
        assert converted.get_quality_string('ascii') == self.solexa_quality
    
    def test_convert_solexa_to_sanger(self):
        """[fastq] convert solexa to sanger"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('solexa')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.solexa_quality)
        converted = self.observed.convert_read_to_format('sanger')
        assert(converted.quality == self.expected_quality).all()
        assert converted.get_quality_string('ascii') == self.sanger_quality
    
    def test_convert_solexa_to_illumina(self):
        """[fastq] convert solexa to illumina"""
        self.observed = fastq.fastqSequencingRead.get_class_by_format('solexa')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.solexa_quality)
        converted = self.observed.convert_read_to_format('illumina')
        assert (converted.quality == self.expected_quality).all()
        assert converted.get_quality_string('ascii') == self.illumina_quality

class TestFastqReader(unittest.TestCase):
    def setUp(self):
        # switch to this directory - so we can have access to data
        os.chdir(os.path.dirname(os.path.abspath( __file__ )))
        seq = 'test-data/sequence.fastq'
        self.seq = fastq.fastqReader(seq)

    def runTest(self):
        """[fastq] reader"""
        for k,v in enumerate(self.seq):
            assert v.identifier == fastq_sequence[k][0]
            assert v.sequence == fastq_sequence[k][1]
            assert v.get_quality_string('ascii') == fastq_sequence[k][2]
            assert (v.quality == fastq_sequence[k][3]).all()

class TestFastqMethods(unittest.TestCase):
    def setUp(self):
        self.s = fastq.fastqSequencingRead()
        self.s.identifier = "@chr5_6255117_6255601_0:0:0_1:0:0_13/1"
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
        self.s.set_quality("&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$")
    
    def test_snapshot(self):
        """[fastq] snapshot the untrimmed sequence"""
        self.s.trim(5)
        self.failIfEqual(self.s.sequence, self.s.sequence_snapshot)
        self.failIfEqual(self.s.quality, self.s.quality_snapshot)
        assert self.s.sequence_snapshot == "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
        assert (self.s.quality_snapshot == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
        assert self.s.get_quality_string('ascii') == "&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&"
    
    def test_trim_5_bases(self):
        """[fastq] trim bases w/ quality < 5"""
        self.s.trim(5)
        assert self.s.sequence == "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGA"
        correct_qual_trim = numpy.array([ 5,  5, 17, 17, 3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 5])
        assert (self.s.quality == correct_qual_trim).all()
        assert (self.s.quality == self.s.quality_snapshot[:-2]).all()
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        assert self.s.trimming == 't'
        
    def test_trim_10_bases(self):
        """[fastq] trim bases w/ quality < 10"""
        self.s.trim(10)
        assert self.s.sequence == "TGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAG"
        correct_qual_trim = numpy.array([ 17, 17, 3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17])
        assert (self.s.quality == correct_qual_trim).all()
        assert (self.s.quality == self.s.quality_snapshot[2:-3]).all()
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        assert self.s.trimming == 't'
    
    def test_trim_no_bases(self):
        """[fastq] don't trim below min_qual"""
        self.s.trim(0)
        assert self.s.sequence == "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
        assert (self.s.quality == numpy.array([ 5,  5, 17, 17, 3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 5, 4, 3])).all()
        assert self.s.get_quality_string('ascii') == "&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$"
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        self.assertFalse(self.s.trimming)
    
    def test_mask_5_bases(self):
        """[fastq] mask bases w/ quality < 5"""
        self.s.mask(5)
        assert self.s.sequence == "CTTGNATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGANN"
        assert (self.s.quality == numpy.array([ 5,  5, 17, 17,  0, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  0,  0])).all()
        assert self.s.sequence.count('N') == 3
        for i in [4, 98, 99]:
            assert self.s.sequence[i] == 'N'
            assert self.s.quality[i] == 0
        assert len(numpy.where(self.s.quality == 0)[0]) == 3
        assert self.s.trimming == 'm'
    
    def test_mask_and_trim_5_bases(self):
        """[fastq] mask and trim bases (edges) where quality < 5"""
        self.s.mask_and_trim(5)
        assert self.s.sequence == "CTTGNATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGA"
        assert (self.s.quality == numpy.array([ 5,  5, 17, 17,  0, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
               17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5])).all()
        assert self.s.quality[4] == 0
        assert len(numpy.where(self.s.quality == 0)[0]) == 1
        assert self.s.trimming == 'mt'
    
    def test_fastq_clone_and_slice(self):
        """[fastq] clone and slice"""
        new = self.s.slice(10, 20)
        assert new.sequence == 'ATGAAAATGC'
        assert (new.quality == numpy.array([17, 17, 17, 17, 17, 17, 17, 17, 17, 17])).all()
        # check to make sure object is cloned
        self.assertNotEqual(id(self.s), id(new))
    
    def test_fastq_no_clone_and_slice(self):
        """[fastq] slice without cloning"""
        new = self.s.slice(10, 20, False)
        assert new.sequence == 'ATGAAAATGC'
        assert (new.quality == numpy.array([17, 17, 17, 17, 17, 17, 17, 17, 17, 17])).all()
        # check to make sure object is cloned
        self.assertEqual(id(self.s), id(new))
    
class TestSangerQualitySpecMethods(unittest.TestCase):
    
    def setUp(self):
        self.s = fastq.fastqSequencingRead()
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
    
    def test_fastq_is_valid_quality_format(self):
        """[fastq] test valid sanger quality format"""
        self.s.set_quality('&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')
        self.failUnless(self.s.quality.any())
        
    def test_fastq_invalid_format_fails(self):
        """[fastq] test invalid quality formats (space in fastq)"""
        self.failUnlessRaises(AssertionError, self.s.set_quality, ' &22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')

class TestNonSangerQualitySpecMethods(unittest.TestCase):    
    
    def setUp(self):
        self.s = fastq.fastqSequencingRead.get_class_by_format('solexa')()
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
    
    def test_fastq_valid_quality_format(self):
        """[fastq] test valid solexa quality format"""
        self.s.set_quality('CCQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQCB@')
        self.failUnless(self.s.quality.any())
    
    def test_fastq_invalid_format_fails2(self):
        """[fastq] test invalid quality formats (incorrect ascii character)"""
        self.failUnlessRaises(AssertionError, self.s.set_quality, '&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')

        
fastq_sequence = (
        ("@chr5_6255117_6255601_0:0:0_1:0:0_13/1", 
            "TTTGACAACTGTTGTTACCCTCCTGTTTATCATGGAGATATCTTTTTTTCTGTATGAATGGCCACAAATGGTCAGACAGAATAGGATGTAGCCTGGAGTA",
            "2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222",
            numpy.array([17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17])
            ),
        ("@chr5_56045273_56045722_1:0:0_2:0:0_12/1",
            "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT",
            "&&22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$",
            numpy.array([ 5,  5, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                   17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])
        )
    )


if __name__ == '__main__':
    unittest.main()
