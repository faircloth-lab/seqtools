#!/usr/bin/env python
# encoding: utf-8

"""
test_fasta.py

Created by Brant Faircloth on 14 December 2010 22:16 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import os
import numpy
import shutil
import unittest
import tempfile
from seqtools.sequence import fasta

#import pdb

class TestSequenceSetQuality(unittest.TestCase):
    def setUp(self):
        self.s = fasta.FastaSequence()
        self.s.identifier = 'chr5_6255117_6255601_0:0:0_1:0:0_13'
        self.s.sequence = 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAATCTGGCAAAGAGCCTACGTGTATTGTGTCCAGTGGGAACAATGCTATGTCACCGAGTCTGTAAGAAT'
    
    def test_set_quality(self):
        """[fasta] convert string quality to array"""
        q = "5 5 17 17 3 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
        +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
        +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
        +" 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17"\
        +" 17 17 17 17 17 17 17 17 17 5 4 3"
        self.s.set_quality(q)
        assert (self.s.quality == numpy.array([ 5,  5, 17, 17,  3, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
                       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,  5,  4,  3])).all()

    def tearDown(self):
        pass

class TestFastaReader(unittest.TestCase):
    def setUp(self):
        # switch to this directory - so we can have access to data
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        seq = 'test-data/sequence.fasta'
        self.fasta = fasta.FastaReader(seq)
    
    def runTest(self):
        """[fasta] reader"""
        for k,v in enumerate(self.fasta):
            assert v.identifier == fasta_sequence[k][0]
            assert v.sequence == fasta_sequence[k][1]
            
    def tearDown(self):
        self.fasta.close()

class TestFastaQualReader(unittest.TestCase):
    def setUp(self):
        # switch to this directory - so we can have access to data
        os.chdir(os.path.dirname(os.path.abspath( __file__ )))
        seq = 'test-data/sequence.fasta'
        qual = 'test-data/sequence.qual'
        self.seq = fasta.FastaQualityReader(seq, qual)

    def runTest(self):
        """[fasta] qual reader"""
        for k,v in enumerate(self.seq):
            assert v.identifier == fasta_sequence[k][0]
            assert v.sequence == fasta_sequence[k][1]
            assert (v.quality == fasta_sequence[k][2]).all()
    
    def tearDown(self):
        self.seq.close()

class TestFastaWriter(unittest.TestCase):
    def setUp(self):
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        seq = 'test-data/sequence.fasta'
        qual = 'test-data/sequence.qual'
        self.seq = fasta.FastaQualityReader(seq, qual)
    
    def _read_raw_contents(self, file):
        return open(file).read()
    
    def test_fasta_write(self):
        """[fasta] fasta writing"""
        d = tempfile.mkdtemp()
        outf = fasta.FastaWriter(os.path.join(d,'test_write.fasta'))
        for s in self.seq:
            outf.write(s)
        outf.close()
        old = self._read_raw_contents('test-data/sequence.fasta')
        new = self._read_raw_contents(os.path.join(d,'test_write.fasta'))
        assert old == new
        shutil.rmtree(d)
    
    def test_fasta_qual_write(self):
        """[fasta] fasta+qual writing"""
        d = tempfile.mkdtemp()
        f = os.path.join(d, 'test_write.fasta')
        q = os.path.join(d, 'test_write.qual')
        outf = fasta.FastaWriter(f,q)
        for s in self.seq:
            outf.write(s)
        outf.close()
        old_s = self._read_raw_contents('test-data/sequence.fasta')
        old_q = self._read_raw_contents('test-data/sequence.qual')
        new_s = self._read_raw_contents(f)
        new_q = self._read_raw_contents(q)
        assert old_s == new_s
        assert old_q == new_q
        shutil.rmtree(d)

fasta_sequence = \
    ((">NoMID rank=0000030 x=1238.5 y=1737.0 length=47","AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG", 
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40])),
    (">MID15_But_Qual_Trimmed rank=0000030 x=1238.5 y=1737.0 length=47","ATACGACGTAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([ 9,  9,  9,  9, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_ATACGACGTA rank=0000030 x=1238.5 y=1737.0 length=47","ATACGACGTAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_OneError_TTACGACGTA rank=0000030 x=1238.5 y=1737.0 length=47","TTACGACGTAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_OneErrorOneGap_TTACGACGT rank=0000030 x=1238.5 y=1737.0 length=47","TTACGACGTGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_NoError rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_OneError rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACCTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_TwoError rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACCTCCTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_OneError_SimpleX1_OneError rank=0000030 x=1238.5 y=1737.0 length=62","ATACCACGTAACCTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_OneGap rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACGTCGTGCGGATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_TwoGap rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACTCGTGCGGATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_NoError_FandR rank=0000030 x=1238.5 y=1737.0 length=77","ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATTCCGCACGACGT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_OneError_Forward_FandR rank=0000030 x=1238.5 y=1737.0 length=77","ATACGACGTAACGTCGTGCGGTATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATTCCGCACGACGT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_OneError_Reverse_FandR rank=0000030 x=1238.5 y=1737.0 length=77","ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATGCCGCACGACGT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_OneError_Forward_OneError_Reverse_FandR rank=0000030 x=1238.5 y=1737.0 length=77","ATACGACGTAACGTCGTGCGGTATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATTGCGCACGACGT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID_NoError_SimpleX_OneErrorOneGap_short_alignment rank=0006231 x=1394.0 y=2171.5 length=51","ATACGACGTAACGTCGTGCGGAGATGTGTATGGGATGTATGTAGGATGTGT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_NoError_F_NEQ_R rank=0000030 x=1238.5 y=1737.0 length=77","ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGATTCCGCTGCTGCT",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])),
    (">MID15_NoError_SimpleX1_NoError_SimpleX1_CONCAT rank=0000030 x=1238.5 y=1737.0 length=62","ATACGACGTAACGTCGTGCGGAATCGAGAGAGAGAGAGAGAGACGTCGTGCGGAATCAGAGAGAGAGAGAGAGAGAG",
        numpy.array([40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
               40, 40, 40, 40, 40, 40, 40, 40, 40])))

if __name__ == '__main__':
    unittest.main()