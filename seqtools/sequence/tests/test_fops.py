#!/usr/bin/env python
# encoding: utf-8

"""
test_fops.py

Created by Brant Faircloth on 14 December 2010 22:15 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import os
import math
import unittest
import multiprocessing
from seqtools.sequence import fops

import pdb

class TestFopsFunctions(unittest.TestCase):
    def test_chunking_by_pieces(self):
        """[fops] chunk by pieces"""
        input = 'test-data/galGal3.read1.fa'
        f_type, delim = fops.file_type(input)
        chunks = fops.get_chunks(input, delim, split_type='pieces')
        tf = fops.make_chunks(chunks, mp=False)
        assert len(tf) == multiprocessing.cpu_count() - 1
        self.clean(tf)
    
    def test_chunking_by_size(self):
        """[fops] chunk by size"""
        input = 'test-data/galGal3.read1.fa'
        f_type, delim = fops.file_type(input)
        chunks = fops.get_chunks(input, delim, mb = 2, split_type='size')
        tf = fops.make_chunks(chunks, mp=False)
        # make sure chunks are ~ mb in size.  problem is that the last one will almost alway be off
        # so just round up and check if w/in 1 mb.
        for f in tf:
            sz = os.path.getsize(f)/1024.**2
            self.failUnlessAlmostEqual(math.ceil(sz), 2, delta=1)
        self.clean(tf)
    
    def test_chunking_by_multiprocessing(self):
        """[fops] chunk by multiprocessing"""
        input = 'test-data/galGal3.read1.fa'
        f_type, delim = fops.file_type(input)
        chunks = fops.get_chunks(input, delim, split_type='pieces')
        tf = fops.make_chunks(chunks, mp=True)
        assert len(tf) == multiprocessing.cpu_count() - 1
        self.clean(tf)
    
    def get_chunks(self, input):
        f_type, delim = fops.file_type(input)
        chunks = fops.get_chunks(input, delim, split_type='pieces')
        values = fops.make_chunks(chunks, mp=False)
        return values
    
    def get_fasta_headers(self, input, headers = set([])):
        for line in open(input, 'rU'):
            if line.startswith('@'):
                headers.add(line.lstrip('@'))
        return headers
    
    def test_fasta_file_type(self):
        """[fops] fasta file type and delimiter check"""
        input = 'test-data/galGal3.read1.fa'
        ft, delim = fops.file_type(input)
        assert [ft, delim] == ['fasta','>']
    
    def test_fastq_file_type(self):
        """[fops] fastq file type and delimiter check"""
        input = 'test-data/galGal3.read1.fq'
        ft, delim = fops.file_type(input)
        assert [ft, delim] == ['fastq','@']
    
    def test_fasta_chunking(self):
        """[fops] sum(chunked fasta file) = unchunked fasta file"""
        input = 'test-data/galGal3.read1.fa'
        unchunk_headers = self.get_fasta_headers(input)
        tf = self.get_chunks(input)
        chunk_headers = set([])
        for f in tf:
            chunk_headers = self.get_fasta_headers(f, chunk_headers)
        assert unchunk_headers == chunk_headers
        self.clean(tf)
    
    def test_fastq_chunking(self):
        """[fops] sum(chunked fastq file) = unchunked fastq file"""
        input = 'test-data/galGal3.read1.fq'
        unchunk_headers = self.get_fasta_headers(input)
        tf = self.get_chunks(input)
        chunk_headers = set([])
        for f in tf:
            chunk_headers = self.get_fasta_headers(f, chunk_headers)
        #pdb.set_trace()
        assert unchunk_headers == chunk_headers
        self.clean(tf)
    
    def clean(self, tf):
        for f in tf:
            os.remove(f)


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    unittest.main()