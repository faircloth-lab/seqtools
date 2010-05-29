#!/usr/bin/env python
# encoding: utf-8

"""
tag_maker.py

Created by Brant Faircloth on 28 May 2010 23:27 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import os
import re
import sys
import string
import itertools
import optparse
import numpy
from levenshtein import getDistance
from levenshtein import getDistanceC


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--tag-length', dest = 'tl', action='store', 
type='int', default = 6, help='The desired tag length')
    
    p.add_option('--edit-distance', dest = 'ed', action='store', 
type='int', default = 3, help='The desired edit distance')

    p.add_option('--no-polybase', dest = 'polybase', action='store_true', default=False, 
help='Remove tags with > 2 identical nucleotides in a row')

    p.add_option('--gc', dest = 'gc', action='store_true', default=False, 
help='Remove tags with GC content (%) 40 > x > 60')

    p.add_option('--comp', dest = 'comp', action='store_true', default=False, 
help='Remove tags that are perfect self-complements')
    
    p.add_option('--use-c', dest = 'clev', action='store_true', default=False, 
help='Use the C version of Levenshtein (faster)')


    (options,arg) = p.parse_args()
    if not options.tl:
        p.print_help()
        sys.exit(2)
    return options, arg

def self_comp(seq):
    '''Return reverse complement of seq'''
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq[::-1].translate(bases)

def main():
    print '''
########################################
#  Tag Generator                       #
########################################
        '''
    count = 0
    good_tags = []
    good_tags_dict = {}
    options, args = interface()
    print '[1] Generating all combinations of tags...'
    all_tags = itertools.product('ACGT', repeat = options.tl)
    if options.polybase:
        regex = re.compile('A{3,}|C{3,}|T{3,}|G{3,}')
    print '[2] If selected, removing tags based on filter criteria...'
    for tag in all_tags:
        #pdb.set_trace()
        tag_seq  = ''.join(tag)
        good = False
        if options.polybase:
            polybase = regex.search(tag_seq)
            if not polybase:
                good = True
        if good and options.gc:
            gc = (tag_seq.count('G') + tag_seq.count('C')) / float(len(tag))
            if 0.40 <= gc <= 0.60:
                good = True
            else:
                good = False
        if good and options.comp:
            if tag_seq != self_comp(tag_seq):
                good = True
            else:
                good = False
        if good:
            tag_name = '{0}'.format(count)
            good_tags.append((tag_name, tag_seq))
            count += 1
    # index the boogers, so we can pull out the good ones
    for k, v in good_tags:
        good_tags_dict.setdefault(k, []).append(v)
    # get the pairwise distance
    print '[3] Calculating the Levenshtein distance across remaining pairs... (Slow)'
    if options.clev:
        print '[C] Using the C version of Levenshtein...'
        distance_pairs = getDistanceC(good_tags, distances = True)
    else:
        distance_pairs = getDistance(good_tags, distances = True)
    hm = numpy.zeros( (count, count) )
    # jam all of our distances into an array
    for tag_pair in distance_pairs:
        t1, t2, d = tag_pair
        hm[t1,t2] = d
    print '[4] Empirically determining the greatest number of returnable tags of Levenshtein distance {0}...'.format(options.ed)
    # get the vector (row) of data with the largest count of tags with edit
    # distance > options.ed
    max_row_sum = [None, 0]
    for pos, row in enumerate(hm):
        bool_row_sum = sum(row >= options.ed)
        if bool_row_sum > max_row_sum[1]:
            max_row_sum = [pos, bool_row_sum]
    # get the indices of those values with edit distances >= options.ed
    best_row = hm[max_row_sum[0]]
    #pdb.set_trace()
    index = numpy.argwhere(best_row >= options.ed)
    keepers = numpy.array([], dtype=int)
    for i in index:
        if not keepers.any():
            keepers = numpy.append(keepers,int(i))
        else:
            # get the column of edit distances for i and reshape it
            column_to_row = numpy.reshape(hm[:,i], len(hm[:,i]))
            # reindex by the keepers
            if sum(column_to_row[keepers] >= options.ed) == len(column_to_row[keepers]):
                keepers = numpy.append(keepers,int(i))
    print '\n'
    for k in keepers:
        print 'Tag{0} = {1}'.format(k, good_tags_dict[str(k)][0]) 
    #pdb.set_trace()    
    #outp.close()

if __name__ == '__main__':
    main()