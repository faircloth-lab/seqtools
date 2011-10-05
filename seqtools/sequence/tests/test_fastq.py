#!/usr/bin/env python
# encoding: utf-8

"""tests for fastq classes and methods

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

import os
import math
import numpy
import cPickle
import unittest
from seqtools.sequence import fastq

#import pdb

class TestSangerQualitySpecMethods(unittest.TestCase):
    """Ensure that invalid Sanger quality values fail"""
    def setUp(self):
        self.s = fastq.FastqSequencingRead()
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
    
    def test_fastq_is_valid_quality_format(self):
        """[fastq] valid sanger quality format passes"""
        self.s.set_quality('&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')
        self.failUnless(self.s.quality.any())
        
    def test_fastq_invalid_format_fails(self):
        """[fastq] invalid quality format (space in fastq) fails"""
        self.failUnlessRaises(AssertionError, self.s.set_quality, ' &22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')

class TestNonSangerQualitySpecMethods(unittest.TestCase):
    """Ensure that invalid Solexa quality values fail"""
    def setUp(self):
        self.s = fastq.FastqSequencingRead.get_class_by_format('solexa')()
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
    
    def test_fastq_valid_quality_format(self):
        """[fastq] valid solexa quality format passes"""
        self.s.set_quality('CCQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQCB@')
        self.failUnless(self.s.quality.any())
    
    def test_fastq_invalid_format_fails2(self):
        """[fastq] invalid solexa quality format (incorrect ascii character) fails"""
        self.failUnlessRaises(AssertionError, self.s.set_quality, '&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$')

class TestFastqClassMethods(unittest.TestCase):
    """Genral test of fastq methods"""
    def setUp(self):
        self.s = fastq.FastqSequencingRead()
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
        """[fastq] ensure that scores do not exceed the maximum for the format"""
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
    
    def test_get_ascii_quality_scores_as_string(self):
        """[fastq] convert decimal scores to ascii scores"""
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
        q = '10 10 20 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 '
        self.s.set_quality(q)
        old = ''.join(self.old_get_ascii_quality_scores(q))
        new = self.s.get_quality_string('ascii')
        assert new == old
    
    def test_get_decimal_quality_scores_as_string(self):
        """[fastq] return decimal quality scores as a string"""
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
        q = '10 10 20 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 '
        self.s.set_quality(q)
        new = self.s.get_quality_string()
        assert new == q.strip()
        
    def test_convert_ascii_to_decimal_quality_scores(self):
        """[fastq] convert ascii scores to decimal scores"""
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
        self.s.set_quality('++5?IIIIIIIIIIIIIII')
        assert (self.s.quality == [10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]).all()
        self.s.set_quality('10 10 20 30 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 ')
        assert (self.s.quality == [10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]).all()
    
    def test_get_sequence(self):
        """[fastq] return sequence"""
        s = self.s.get_sequence()
        assert s == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
    
    def test_is_valid_quality_length(self):
        """[fastq] quality length matches sequence length"""
        self.assertTrue(self.s.is_valid_quality_length())
    
    def  test_is_not_valid_quality_length(self):
        """[fastq] invalid quality length fails"""
        self.s.quality = numpy.array([ 5,  5, 17, 17, 17, 17])
        self.assertFalse(self.s.is_valid_quality_length())
    
    def  test_assert_sequence_quality_lengths(self):
        """[fastq] quality length matches sequence length"""
        self.s.quality = numpy.array([ 5,  5, 17, 17, 17, 17])
        self.failUnlessRaises(AssertionError, self.s.assert_sequence_quality_lengths)
        
class TestFastqErrors(unittest.TestCase):
    def setUp(self):
        self.s = fastq.FastqSequencingRead()
        self.s.sequence = 'CTTGGATCAGATGAAAATG'
    
    def test_get_decimal_quality_when_no_quality_array(self):
        """[fastq] raise ValueError when no quality string"""
        self.assertRaises(ValueError, self.s.get_quality_string)

class TestFastqConversions(unittest.TestCase):
    """test conversions from one format (sanger) to another (solexa)"""
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
        self.observed = fastq.FastqSequencingRead.get_class_by_format('illumina')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.illumina_quality)
        self.assertRaises(AssertionError, self.observed.convert_read_to_format, 'bob')

    def test_convert_ilumina_to_sanger(self):
        """[fastq] convert illumina to sanger"""
        self.observed = fastq.FastqSequencingRead.get_class_by_format('illumina')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.illumina_quality)
        # the array value of this should be equal to our "usual" array value
        assert (self.observed.quality == self.expected_quality).all()
        converted = self.observed.convert_read_to_format('sanger')
        assert converted.get_quality_string('ascii') == self.sanger_quality
    
    def test_convert_sanger_to_illumina(self):
        """[fastq] convert sanger to illumina"""
        self.observed = fastq.FastqSequencingRead.get_class_by_format('sanger')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.sanger_quality)
        assert (self.observed.quality == self.expected_quality).all()
        converted = self.observed.convert_read_to_format('illumina')
        assert converted.get_quality_string('ascii') == self.illumina_quality
    
    def test_convert_sanger_to_solexa(self):
        """[fastq] convert sanger to solexa"""
        self.observed = fastq.FastqSequencingRead.get_class_by_format('sanger')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.sanger_quality)
        converted = self.observed.convert_read_to_format('solexa')
        assert (converted.quality == self.expected_solexa_quality).all()
        assert converted.get_quality_string('ascii') == self.solexa_quality
    
    def test_convert_solexa_to_sanger(self):
        """[fastq] convert solexa to sanger"""
        self.observed = fastq.FastqSequencingRead.get_class_by_format('solexa')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.solexa_quality)
        converted = self.observed.convert_read_to_format('sanger')
        assert(converted.quality == self.expected_quality).all()
        assert converted.get_quality_string('ascii') == self.sanger_quality
    
    def test_convert_solexa_to_illumina(self):
        """[fastq] convert solexa to illumina"""
        self.observed = fastq.FastqSequencingRead.get_class_by_format('solexa')()
        self.observed.identifier = self.identifier
        self.observed.sequence = self.sequence
        self.observed.set_quality(self.solexa_quality)
        converted = self.observed.convert_read_to_format('illumina')
        assert (converted.quality == self.expected_quality).all()
        assert converted.get_quality_string('ascii') == self.illumina_quality

class TestFastqSummarizer(unittest.TestCase):
    """test our methods to summarize fastq read data"""
    def setUp(self):
        fastqs = fastq.FastqReader('../tests/test-data/sequence_for_summarizer.fastq')
        self.summary = fastq.FastqSummarizer()
        for read in fastqs:
            self.summary.add(read)
        self.archive = None
        
    def tearDown(self):
        if self.archive:
            self.archive.close()
    
    def test_correct_quality_stack(self):
        """[fastq] summary stack qualities"""
        archive = open('../tests/test-data/sequence_for_summarizer_quality.pickle','r')
        archive_data = cPickle.loads(archive.read())
        # this is a little kludgy - but drop NaNs and compare.  The main
        # problem is that nan values are not equiv.
        # >>> numpy.nan == numpy.nan
        # False
        assert (archive_data[numpy.isfinite(archive_data)] \
            == self.summary.qualities[numpy.isfinite(self.summary.qualities)]).all()
        
    def test_correct_base_counts(self):
        """[fastq] summary count bases"""
        archive = open('../tests/test-data/sequence_for_summarizer_bases.pickle','r')
        archive_data = cPickle.loads(archive.read())
        assert (archive_data == self.summary.bases).all()
    
    def test_correct_lengths(self):
        """[fastq] summary sequence lengths"""
        expected = {100:2, 200:1}
        assert self.summary.lengths == expected
    
    def test_correct_decimal_range(self):
        """[fastq] summary decimal range"""
        assert self.summary.get_decimal_range() == (3.0, 17.0)
    
    def test_correct_length_count(self):
        """[fastq] summary sequence lengths"""
        expected = {100:2, 200:1}
        assert self.summary.get_length_counts() == expected
    
    def test_max_read_length(self):
        """[fastq] summary max read length"""
        assert self.summary.get_max_read_length() == 200
    
    def test_min_read_lenght(self):
        """[fastq] summary min read length"""
        assert self.summary.get_min_read_length() == 100
    
    def test_read_count_for_column_three(self):
        """[fastq] summary read count for column 10"""
        column_ten = self.summary.get_read_count_for_column(10)
        assert column_ten == 3
    
    def test_read_count_for_column_one_o_one(self):
        """[fastq] summary read count for column 101"""
        column_one_o_one = self.summary.get_read_count_for_column(101)
        assert column_one_o_one == 1
    
    def test_read_count_for_column_two_o_one(self):
        """[fastq] summary read count for column 201"""
        column_too_many = self.summary.get_read_count_for_column(201)
        assert column_too_many == 0
    
    def test_read_count(self):
        """[fastq] summary read count (overall)"""
        assert self.summary.get_read_count() == 3
    
    def test_base_counts_for_column(self):
        """[fastq] summary base counts for column"""
        expected = {'A': 0.0, 'C': 0.0, 'T': 3.0, 'G': 0.0, 'N': 0.0}
        assert self.summary.get_base_counts_for_column(1) == expected
    
    def test_all_scores_for_column(self):
        """[fastq] summary all scores for column"""
        expected = numpy.array([ 17.,   5.,  17.], dtype=numpy.float32)
        assert (self.summary.get_all_scores_for_column(1) == expected).all()
    
    def test_finite_scores_for_column(self):
        """[fastq] summary finite scores for column"""
        expected = numpy.array([ 17.], dtype=numpy.float32)
        assert (self.summary.get_finite_scores_for_column(101) == expected).all()

    def test_nan_count_for_column(self):
        """[fastq] summary NaN count for column"""
        assert self.summary.get_nan_count_for_column(101) == 2
        
    def test_quality_min_for_column(self):
        """[fastq] summary quality min for column"""
        assert self.summary.get_quality_min_for_column(1) == 5.0

    def test_quality_max_for_column(self):
        """[fastq] summary quality max for column"""
        assert self.summary.get_quality_max_for_column(1) == 17.0
    
    def test_quality_average_for_column(self):
        """[fastq] summary average quality for column"""
        assert self.summary.get_quality_average_for_column(1) == 13.0
    
    def test_quality_std_deviation_for_column(self):
        """[fastq] summary quality standard deviation for column"""
        self.assertAlmostEqual(self.summary.get_quality_std_deviation_for_column(1), 6.9282032302755088, 2)
    
    def test_summary_for_column(self):
        """[fastq] summary for column"""
        expected = {'quality_stddev': 6.9282032302755088, 'column': 1, 'quality_avg': 13.0, 'quality_max': 17.0, 'quality_min': 5.0, 'reads': 3, 'nans': 0}
        observed = self.summary.get_summary_for_column(1)
        for element in observed:
            self.assertAlmostEqual(observed[element], expected[element], 2)
        
class TestFastqReader(unittest.TestCase):
    """test ability to parse fastq files"""
    def setUp(self):
        # switch to this directory - so we can have access to data
        os.chdir(os.path.dirname(os.path.abspath( __file__ )))
        seq = 'test-data/sequence.fastq'
        self.seq = fastq.FastqReader(seq)

    def runTest(self):
        """[fastq] reader"""
        for k,v in enumerate(self.seq):
            assert v.identifier == fastq_sequence[k][0]
            assert v.sequence == fastq_sequence[k][1]
            assert v.get_quality_string('ascii') == fastq_sequence[k][2]
            assert (v.quality == fastq_sequence[k][3]).all()

class TestFastqWriter(unittest.TestCase):
    """test ability to write fastq files"""
    def setUp(self):
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        fastq_file = 'test-data/sequence.fastq'
        self.seq = fastq.FastqReader(fastq_file)
    
    def _read_raw_contents(self, file):
        return open(file).read()
    
    def test_fastq_write_sanger_format(self):
        """[fastq] fastq writing (sanger format)"""
        outf = fastq.FastqWriter('test-output/test_write.fastq')
        for s in self.seq:
            outf.write(s)
        outf.close()
        old = self._read_raw_contents('test-data/sequence.fastq')
        new = self._read_raw_contents('test-output/test_write.fastq')
        assert old == new
    
    def test_fastq_write_solexa_format(self):
        """[fastq] fastq writing (solexa format)"""
        outf = fastq.FastqWriter('test-output/test_write.fastq', format='solexa')
        for s in self.seq:
            outf.write(s)
        outf.close()
        old = self._read_raw_contents('test-data/sequence.solexa.fastq')
        new = self._read_raw_contents('test-output/test_write.fastq')
        assert old == new
    
    def test_fastq_write_illumina_format(self):
        """[fastq] fastq writing (illumina format)"""
        outf = fastq.FastqWriter('test-output/test_write.fastq', format='illumina')
        for s in self.seq:
            outf.write(s)
        outf.close()
        old = self._read_raw_contents('test-data/sequence.illumina.fastq')
        new = self._read_raw_contents('test-output/test_write.fastq')
        assert old == new

    def tearDown(self):
        os.remove('test-output/test_write.fastq')

class TestFastqTrimmingMaskingAndSlicingMethods(unittest.TestCase):
    """test snapshotting, trimming, masking, masking and trimmming, and slicing"""
    def setUp(self):
        self.s = fastq.FastqSequencingRead()
        self.s.identifier = "@chr5_6255117_6255601_0:0:0_1:0:0_13/1"
        self.s.sequence = "CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT"
        self.s.set_quality("&&22$22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222&%$")
    
    def test_snapshot(self):
        """[fastq] snapshot the untrimmed sequence"""
        self.s.trim(5, clone = False)
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
        self.s.trim(5, clone = False)
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
        self.s.trim(10, clone = False)
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
        self.s.trim(0, clone = False)
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
        self.s.mask(5, clone = False)
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
        self.s.mask_and_trim(5, clone = False)
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
