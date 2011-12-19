from math import log
from collections import namedtuple

import pdb

class PslReader():
    """Represents an iterator over psl files from blat"""

    def __init__(self, psl):
        """set the psl file attribute"""
        self.psl = open(psl)
        self.psl_template = namedtuple('Psl',
            [
                'match',
                'mismatch',
                'repmatch',
                'N',
                'q_gap_count',
                'q_gap_bases',
                't_gap_count',
                't_gap_bases',
                'strand',
                'q_name',
                'q_size',
                'q_start',
                'q_end',
                't_name',
                't_size',
                't_start',
                't_end',
                'block_count',
                'block_sizes',
                'qstarts',
                't_starts',
            ])

    def close(self):
        """close file"""
        self.psl.close()

    def next(self):
        line = self.psl.readline()
        if not line:
            raise StopIteration
        # strip and split
        line = line.strip().split('\t')
        line = [int(val) if val.isdigit() else val for k,val in enumerate(line)]
        return self.psl_template._make(line)

    def __iter__(self):
        """iterator"""
        while True:
            yield self.next()


def percent_id(psl, protein=False, mRNA=True):
    '''
    convert the percent sequence ID of an alignment from a single line of a 
    parsed PSL file.  Code adapted from 
    http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    millibad = 0
    if protein:
        sizeMul = 3
    else:
        sizeMul = 1
    qAliSize = sizeMul * (psl.q_end - psl.q_start)
    tAliSize = psl.t_end - psl.t_start
    aliSize = min(qAliSize, tAliSize)
    if aliSize <= 0:return 0
    else:
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
            if mRNA:
                sizeDif = 0;
            else:
                sizeDif = -sizeDif
        insertFactor = psl.q_gap_count
        if not mRNA:
            insertFactor += psl.t_gap_count
        total = (sizeMul * (psl.match + psl.repmatch + psl.mismatch))
        if total != 0:
            milliBad = (1000 * (psl.mismatch * sizeMul + insertFactor + round(3*log(1+sizeDif)))) / total
        percent = round(100 - milliBad * 0.1,0)
        return percent

def score(psl, sizeMul = 1):
    '''
    convert the score of an alignment from a single line of a parsed PSL file.
    Code adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4
    '''
    score = sizeMul * (psl.match + (psl.repmatch >> 1)) - sizeMul * psl.mismatch - psl.q_gap_count - psl.t_gap_count
    return score
