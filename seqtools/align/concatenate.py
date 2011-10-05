from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment

class ConcatenatedAlignment(Alignment):
    '''A child class of the biopython alignment class, created so that we could
    add the concat method to concatenate multiple alignments'''
    def __init__(self):
        Alignment.__init__(self, Gapped(IUPAC.unambiguous_dna, '-'))
    
    def add(self, alignment):
        '''concatenate alignment objects (on a per SeqRecord basis)'''
        #pdb.set_trace()
        if not self._records:
            for seq in alignment:
                seq.name = seq.id
                self._records.append(seq)
        else:
            for seq in alignment:
                for pos, record in enumerate(self._records):
                    # this assumes that we will be using id as
                    # the joining attribute...
                    if seq.id == record.id:
                        c_seq       = SeqRecord(record.seq + seq.seq)
                        c_seq.name  = record.id
                        c_seq.id    = record.id
                        self._records[pos] = c_seq