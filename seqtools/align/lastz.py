#!/usr/bin/env python
# encoding: utf-8

"""
Lastz.py

Created by Brant Faircloth on 01 May 2010 19:01 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  This is a helper class to simplify calls to lastz.

USAGE:  

    import Lastz
    
    lastz = Lastz.Align(target, query, matchcount, identity, output)
    lastz_stderr, lastz_stdout = lastz.run()
    results_file = lastz.output

"""

import os
import tempfile
import subprocess
from collections import namedtuple
from collections import defaultdict

#import pdb

class Align():
    '''docstring for lastz'''
    def __init__(self, target, query, matchcount, identity, out=False):
        # if not an output file, create a temp file to hold output
        if not out:
            fd, self.output = tempfile.mkstemp(suffix='.lastz')
            os.close(fd)
        else:
            self.output = out
        self.cli = 'lastz {0}[multiple,nameparse=full] {1}[nameparse=full]\
            --strand=both \
            --seed=12of19 \
            --transition \
            --nogfextend \
            --nochain \
            --gap=400,30 \
            --xdrop=910 \
            --ydrop=8370 \
            --hspthresh=3000 \
            --gappedthresh=3000 \
            --noentropy \
            --coverage={2} \
            --identity={3} \
            --output={4} \
            --format=general-:score,name1,strand1,zstart1,end1,length1,name2,\
strand2,zstart2,end2,length2,diff,cigar,identity,\
continuity'.format(target, query, matchcount, identity, self.output)
    
    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(self.cli, shell=True, \
            stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        return lastz_stdout, lastz_stderr

class Reader():
    """read a lastz file and return an iterator over that file"""
    def __init__(self, lastz_file, long_format = False):
        self.file = open(lastz_file, 'rU')
        self.long_format = long_format
        
    def __del__(self):
        """close files"""
        self.file.close()
        
    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()
        
    def next(self):
        """read next fastq sequence and return as named tuple"""
        lastz_result = self.file.readline()
        if not lastz_result:
            raise StopIteration
        if not self.long_format:
            Lastz = namedtuple('Lastz', 'score,name1,strand1,zstart1,end1,length1,name2,'+
            'strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,'+
            'continuity,percent_continuity')
        else:
            Lastz = namedtuple('Lastz', 'score,name1,strand1,zstart1,end1,length1,name2,'+
            'strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,'+
            'continuity,percent_continuity,coverage,percent_coverage')
        aligns = defaultdict(lambda: defaultdict(list))
        lastz_result_split = lastz_result.strip('\n').split('\t')
        for k,v in enumerate(lastz_result_split):
            if k in [3,4,5,8,9,10]:
                lastz_result_split[k] = int(v)
            elif '%' in v:
                lastz_result_split[k] = float(v.strip('%'))
        return Lastz._make(lastz_result_split)

def get_name(header, splitchar = "_", items = 2):
    """use own function vs. import from match_contigs_to_probes - we don't want lowercase"""
    if splitchar:
        return "_".join(header.split(splitchar)[:items]).lstrip('>')
    else:
        return header.lstrip('>')

def get_query_dupes(lastz):
    assert isinstance(lastz, Reader), "Object must be of lastz.Reader class"
    matches = defaultdict(list)
    dupes = set()
    #import pdb
    #pdb.set_trace()
    for lz in lastz:
        query_name = get_name(lz.name2, "|", 1)
        matches[query_name].append(1)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k,v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            dupes.add(k)
    return dupes

if __name__ == '__main__':
    pass
