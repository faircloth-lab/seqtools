#!/usr/bin/env python
# encoding: utf-8

"""
test_sequence.py

Created by Brant Faircloth on 15 December 2010 13:47 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import os
import numpy
import unittest
from demuxipy.lib import sequence

import pdb

class TestSequenceSetQuality(unittest.TestCase):
    def setUp(self):
        self.s = sequence.SequencingRead()
        self.s.identifier = '@chr5_6255117_6255601_0:0:0_1:0:0_13'
        self.s.sequence = 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
    
    def test_set_quality(self):
        """[sequence] convert string quality to array"""
        q = "5 5 17 17 3 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 5 4 3"
        self.s.set_quality(q)
        assert (self.s.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()

    def tearDown(self):
        pass

class TestSequenceMethods(unittest.TestCase):
        
    def setUp(self):
        self.s = sequence.SequencingRead()
        self.s.identifier = '@chr5_6255117_6255601_0:0:0_1:0:0_13'
        self.s.sequence = 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        q = "5 5 17 17 3 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
            +" 17 17 17 17 17 17 17 17 17 17 17 5 4 3"
        self.s.set_quality(q)
    
    def test_append_sequence(self):
        """[sequence] append sequence"""
        addition = self.s.sequence + 'CTTGGATCAGATGAAAA'
        self.s.append_sequence('CTTGGATCAGATGAAAA')
        assert self.s.sequence == addition
    
    def test_append_quality(self):
        """[sequence] append quality"""
        additions = numpy.concatenate((self.s.quality, numpy.array([10, 10, 10, 10, 10, 20])), axis = 0)
        self.s.append_quality('10 10 10 10 10 20 ')
        assert (self.s.quality == additions).all()
    
    def test_is_dna(self):
        """[sequence] is DNA"""
        assert self.s.is_DNA()
        self.s.sequence = 'CUUGGAUCAGATGAAAAU'
        self.assertFalse(self.s.is_DNA())
    
    def test_is_rna(self):
        """[sequence] is RNA"""
        self.assertFalse(self.s.is_RNA())
        self.s.sequence = 'CUUGGAUCAGATGAAAAU'
        assert self.s.is_RNA()
    
    def test_clone(self):
        """[sequence] clone sequence object"""
        t = self.s.clone()
        self.assertNotEqual(id(self.s), id(t))
    
    def test_snapshot(self):
        """[sequence] snapshot sequence and quality"""
        self.s.snapshot()
        assert self.s.sequence_snapshot == self.s.sequence
        self.s.sequence = None
        assert self.s.sequence_snapshot == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTA'\
            +'TTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        assert (self.s.quality_snapshot == self.s.quality).all()
        self.s.quality = None
        assert (self.s.quality_snapshot == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
    
    def test_reverse(self):
        """[sequence] reverse strand"""
        t = self.s.reverse()
        self.assertNotEqual(id(self.s), id(t))
        assert t.sequence == 'TAAGAATGTCTGAGCCACTGTATCGTAACAAGGGTGACCTGTGTTATGTGCATCCGAGAAACGGTCTA'\
            +'ATTTATGTTCGACGTAAAAGTAGACTAGGTTC'
        assert (t.quality == numpy.array([ 3,  4,  5, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  3, 17, 17,  5,  5])).all()
        # ensure this works without cloning, too
        r = t.reverse(False)
        assert id(t) == id(r)
        assert r.sequence == t.sequence
        assert (r.quality == t.quality).all()
        assert t.sequence == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGG'\
            +'GAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        assert (t.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
    
    def test_DNA_complement(self):
        """[sequence] complement strand"""
        t = self.s.complement()
        self.assertNotEqual(id(self.s), id(t))
        assert t.sequence == 'GAACCTAGTCTACTTTTACGTCGAACATAAATTAGACCGTTTCTCGGATGCACATAACACAGGTCACC'\
            +'CTTGTTACGATACAGTGGCTCAGACATTCTTA'
        assert (t.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
        # ensure this works without cloning, too
        r = t.complement(False)
        assert id(t) == id(r)
        assert r.sequence == t.sequence
        assert (r.quality == t.quality).all()
        assert t.sequence == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGG'\
            +'GAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        assert (t.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
    
    def test_RNA_complement(self):
        """[sequence] complement strand"""
        t = self.s.sequence_as_RNA()
        self.assertNotEqual(id(self.s), id(t))
        t.complement()
        assert t.sequence == 'CUUGGAUCAGAUGAAAAUGCAGCUUGUAUUUAAUCUGGCAAAGAGCCUACGUGUAUUGUGUCCAGUGGG'\
            +'AACAAUGCUAUGUCACCGAGUCUGUAAGAAU'
        assert (t.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
    
    def test_reverse_complement(self):
        """[sequence] reverse complement strand"""
        t = self.s.reverse_complement()
        self.assertNotEqual(id(self.s), id(t))
        assert t.sequence == 'ATTCTTACAGACTCGGTGACATAGCATTGTTCCCACTGGACACAATACACGTAGGCTCTTTGCCAGAT'\
            +'TAAATACAAGCTGCATTTTCATCTGATCCAAG'
        assert (t.quality == numpy.array([ 3,  4,  5, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  3, 17, 17,  5,  5])).all()
        # ensure this works without cloning, too
        r = t.reverse_complement(False)
        assert id(t) == id(r)
        assert r.sequence == t.sequence
        assert (r.quality == t.quality).all()
        assert t.sequence == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGG'\
            +'GAACAATGCTATGTCACCGAGTCTGTAAGAAT'
        assert (t.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                           17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()
    
    def test_DNA_to_RNA_to_DNA(self):
        """[sequence] convert from DNA to RNA to DNA"""
        t = self.s.sequence_as_RNA()
        self.assertNotEqual(id(self.s), id(t))
        assert t.sequence == 'CUUGGAUCAGAUGAAAAUGCAGCUUGUAUUUAAUCUGGCAAAGAGCCUACGUGUAUUGUGUCCAGUGG'\
            +'GAACAAUGCUAUGUCACCGAGUCUGUAAGAAU'
        r = t.sequence_as_DNA(False)
        assert id(t) == id(r)
        assert t.sequence == 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGG'\
            +'GAACAATGCTATGTCACCGAGTCTGTAAGAAT'
            
    
class TestSequenceTrimmingMethods(unittest.TestCase):
    def setUp(self):
        self.s = sequence.SequencingRead()
        self.s.sequence = "ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACGACGT"
        self.s.quality = numpy.array([10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30, 20, 20, 10, 10])

    def test_snapshot(self):
       """[sequence] snapshot untrimmed sequence"""
       self.s.trim(20)
       self.failIfEqual(self.s.sequence, self.s.sequence_snapshot)
       self.failIfEqual(self.s.quality, self.s.quality_snapshot)
       assert self.s.sequence_snapshot == "ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACGACGT"
       assert (self.s.quality_snapshot == numpy.array([10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                  40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                  40, 40, 40, 40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40,
                  40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                  40, 40, 40, 40, 30, 20, 20, 10, 10])).all()

    def test_trim_20_bases(self):
        """[sequence] trim bases w/ quality < 20"""
        self.s.trim(20)
        assert self.s.sequence == "ACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACGAC"
        correct_qual_trim = numpy.array([20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30, 20, 20])
        assert (self.s.quality == correct_qual_trim).all()
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        assert self.s.trimming == 't'

    def test_trim_30_bases(self):
        """[sequence] trim bases w/ quality < 30"""
        self.s.trim(30)
        assert self.s.sequence == "CGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACG"
        correct_qual_trim = numpy.array([30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30])
        assert (self.s.quality == correct_qual_trim).all()
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        assert self.s.trimming == 't'

    def test_trim_no_bases(self):
        """[sequence] trim no bases when the min_qual = 0"""
        self.s.trim(0)
        assert self.s.sequence == "ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACGACGT"
        assert (self.s.quality == numpy.array([10, 10, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 10, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30, 20, 20, 10, 10])).all()
        # make sure we didn't mask anything inside the edges
        assert 'N' not in self.s.sequence
        assert not self.s.trimming

    def test_mask_20_bases(self):
        """[sequence] mask bases w/ quality < 20"""
        self.s.mask(20)
        assert self.s.sequence == "NNACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGANAGAGAGAGAGAGAGAGAGGATGCCGCACGACNN"
        assert (self.s.quality == numpy.array([0, 0, 20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 0, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30, 20, 20, 0, 0])).all()
        for i in [0,1,43,75,76]:
            assert self.s.sequence[i] == 'N'
            assert self.s.quality[i] == 0
        assert self.s.sequence.count('N') == 5
        assert len(numpy.where(self.s.quality == 0)[0]) == 5
        assert self.s.trimming == 'm'

    def test_mask_and_trim_20_bases(self):
        """[sequence] mask and trim bases (edges) where quality < 5"""
        self.s.mask_and_trim(20)
        assert self.s.sequence == "ACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGANAGAGAGAGAGAGAGAGAGGATGCCGCACGAC"
        assert (self.s.quality == numpy.array([20, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 0, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                   40, 40, 40, 40, 30, 20, 20])).all()
        # make sure we've masked the correct base and no other bases
        assert self.s.sequence[41] == 'N'
        assert self.s.sequence.count('N') == 1
        assert self.s.quality[41] == 0
        assert len(numpy.where(self.s.quality == 0)[0]) == 1
        assert self.s.trimming == 'mt'
    
    def test_clone_and_slice(self):
        """[sequence] clone and slice"""
        new = self.s.slice(10, 20)
        assert new.sequence == 'ACGTCGTGCG'
        assert (new.quality == numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40])).all()
        # check to make sure object is cloned
        self.assertNotEqual(id(self.s), id(new))
    
    def test_no_clone_and_slice(self):
        """[sequence] slice without cloning"""
        new = self.s.slice(10, 20, False)
        assert new.sequence == 'ACGTCGTGCG'
        assert (new.quality == numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40])).all()
        # check to make sure object is cloned
        self.assertEqual(id(self.s), id(new))
    
    def tearDown(self):
        pass
    

if __name__ == '__main__':
    unittest.main()