#!/usr/bin/env python
# encoding: utf-8

"""
test_transform.py

Created by Brant Faircloth on 15 December 2010 22:21 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import os
import numpy
import math
import unittest
import exceptions
from demuxipy.lib import transform

import pdb

class TestFastqClassMethods(unittest.TestCase):
    def setUp(self):
        self.sequence = 'CTTGGATCAGATGAAAATGCAGCTTGTATTTAAT'
    
    def test_reverse(self):
        """[transform] reverse"""
        assert transform.reverse(self.sequence) == 'TAATTTATGTTCGACGTAAAAGTAGACTAGGTTC'

    def test_DNA_complement(self):
        """[transform] DNA complement"""
        assert transform.DNA_complement(self.sequence) == 'GAACCTAGTCTACTTTTACGTCGAACATAAATTA'
    
    def test_RNA_complement(self):
        """[transform] RNA complement"""
        assert transform.RNA_complement(self.sequence) == 'GTTCCUTGUCUTCUUUUTCGUCGTTCTUTTTUUT'
    
    def test_DNA_reverse_complement(self):
        """[transform] DNA reverse complement"""
        assert transform.DNA_reverse_complement(self.sequence) == 'ATTAAATACAAGCTGCATTTTCATCTGATCCAAG'

    def test_RNA_reverse_complement(self):
        """[transform] RNA reverse complement"""
        assert transform.RNA_reverse_complement(self.sequence) == 'TUUTTTUTCTTGCUGCTUUUUCTUCUGTUCCTTG'
    
    def test_to_DNA(self):
        """[transform] to DNA"""
        assert transform.to_DNA('CUUGGAUCAGAUGAAAAUGCAGCUUGUAUUUAAU') == self.sequence
    
    def test_to_RNA(self):
        """[transform] to RNA"""
        assert transform.to_RNA(self.sequence) == 'CUUGGAUCAGAUGAAAAUGCAGCUUGUAUUUAAU'
        