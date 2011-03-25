#!/usr/bin/env python
# encoding: utf-8

"""functions to convert DNA sequences

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

import string

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = string.maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = string.maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = string.maketrans( "Tt", "Uu" )
RNA_TO_DNA = string.maketrans( "Uu", "Tt" )

def reverse(sequence):
    """return reverse of the sequence strand as string"""
    return sequence[::-1]
    
def DNA_complement(sequence):
    """return complement of the DNA strand as string"""
    return sequence.translate(DNA_COMPLEMENT)
    
def RNA_complement(sequence):
    """return complement of the RNA strand as string"""
    return sequence.translate(RNA_COMPLEMENT)
    
def DNA_reverse_complement(sequence):
    """return reverse complement of the DNA strand as string"""
    sequence = reverse(sequence)
    return DNA_complement(sequence)
    
def RNA_reverse_complement(sequence):
    """return the reverse complement of the RNA strand as string"""
    sequence = reverse(sequence)
    return RNA_complement( sequence )
    
def to_DNA(sequence):
    """convert RNA strand to DNA, returning string"""
    if 'u' in sequence.lower():
        return sequence.translate( RNA_TO_DNA )
    else: # pragma: no cover
        return sequence
    
def to_RNA(sequence):
    """convert DNA strand to RNA, returning string"""
    if not 'u' in sequence.lower():
        return sequence.translate( DNA_TO_RNA )
    else: # pragma: no cover
        return sequence
