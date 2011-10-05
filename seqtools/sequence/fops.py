#!/usr/bin/env python
# encoding: utf-8

"""
io.py

Created by Brant Faircloth on 11 December 2010 11:28 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

Read and return sequence of various formats.

"""

import os
import tempfile
import multiprocessing

#import pdb

def file_type(input):
    """given an input file, determine the type and return both type and record delimiter (> or @)"""
    name, extension = os.path.splitext(os.path.basename(input))
    fastas = set(['.fsa','.fasta','.fa'])
    fastqs = set(['.fastq','.fq'])
    gffs = set(['.gff'])
    if extension in fastas:
        ft = 'fasta'
        delim = '>'
    elif extension in fastqs:
        ft = 'fastq'
        delim = '@'
    # TODO:  sff ???
    #elif extension in gffs:
    #    ft = 'sff'
    #    delim = None
    else: # pragma: no cover
        raise IOError, "Input file not of correct extension"
    return ft, delim
    
def _get_file_chunks(input, delim, size):
    """given input, record delimiter, and chunk size, yield an iterator contains file seek (start)
    and file read (stop) positions.  Return final position as (6365605, None)."""
    f = open(input)
    while 1:
        start = f.tell()
        f.seek(size, 1)
        line = f.readline()
        if not line:
            break
        # if this isn't a fasta header line, read forward until
        # we get to one
        while not line.startswith(delim):
            line = f.readline()
        else:
            # now that we got to a fasta header, we're at the end.
            # back up the length of the fasta header.
            f.seek(-len(line), 1)
            # tuple up
            yield start, f.tell() - start, input
    # make sure we catch the (start, distance) for the end of the file, too
    yield start, None, input
    f.close()
    
def get_chunks(input, delim, split_type, mb=1, splits = None):
    """return a tuple of file seek (start, distance) positions covering chunks of a file"""
    if split_type == 'size':
        size = mb * (1024**2)
    if split_type == 'pieces':
        if not splits:
            splits = multiprocessing.cpu_count() - 1
        size = int(round((os.path.getsize(input)/float(splits)), 0))
    return _get_file_chunks(input, delim, size)

def _split_file(chunk):
    """function to split a file into appropriate pieces given (start, stop) file seek coords"""
    f = open(chunk[2])
    f.seek(chunk[0])
    if chunk[1]:
        d = f.read(chunk[1])
    else:
        d = f.read()
    td, tf = tempfile.mkstemp(suffix='.splt')
    os.close(td)
    otf = open(tf, 'w')
    otf.write(d)
    otf.close()
    f.close()
    return tf
    
def make_chunks(chunks, pool = None, mp = True):
    """return a list of tempfiles that are the chunked input file"""
    if mp and not pool:
        # create a multiprocessing pool
        procs = multiprocessing.cpu_count() - 1
        pool = multiprocessing.Pool(procs)
        chunks = pool.map(_split_file, chunks)
        # close up the pool if we no longer want to swim
        pool.close()
        pool.join()
    else:
        chunks = map(_split_file, chunks)
    return chunks