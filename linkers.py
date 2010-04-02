#!/usr/bin/env python
# encoding: utf-8
"""
linkers.py

Copyright (c) 2009-2010 Brant C. Faircloth.  All rights reserved.

Provides:
    - parsing and error correction of hierarchically tagged
        next generation sequence reads

Part of mc454:
    - parsing hierarchically tagged sequence reads
    - microsatellite identification
    - sequence pooling/clustering
    - microsatellite primer design

"""

import os, sys, re, pdb, time, numpy, string, MySQLdb, ConfigParser, multiprocessing, cPickle, optparse, progress
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.SeqIO import QualityIO
from Bio.Alphabet import SingleLetterAlphabet

def motd():
    '''Startup info'''
    motd = '''
    ##############################################################
    #                     linkers.py                             #
    # Provides:                                                  #
    #   - parsing and error correction of hierarchically tagged  #
    #     next generation sequence reads                         #
    #                                                            #
    # Part of mc454:                                             #
    #   - parsing hierarchically tagged sequence reads           #
    #   - microsatellite identification                          #
    #   - sequence pooling/clustering                            #
    #   - microsatellite primer design                           #
    #                                                            #
    # Copyright (c) 2009-2010 Brant C. Faircloth                 #
    # 621 Charles E. Young Drive                                 #
    # University of California, Los Angeles, 90095, USA          #
    ##############################################################\n
    '''
    print motd

def revComp(seq):
    '''Return reverse complement of seq'''
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq.translate(bases)[::-1]

def revCompTags(tags):
    '''Return the reverse complements of a tag dictionary'''
    revTags = {}
    for tag in tags:
        revTags[revComp(tag)] = tags[tag]
    return revTags

def trim(sequence, left=None, right=None):
    '''Trim a given sequence given left and right offsets'''
    if left and right:
        sequence = sequence[left:right]
    elif left:
        sequence = sequence[left:]
    elif right:
        sequence = sequence[:right]
    return sequence

def matches(tag, seq_match_span, tag_match_span, allowed_errors):
    '''Determine the gap/error counts for a particular match'''
    # deal with case where tag match might be perfect, but extremely gappy, 
    # e.g. ACGTCGTGCGGA-------------------------ATC
    if tag_match_span.count('-') > allowed_errors or \
        seq_match_span.count('-') > allowed_errors:
        return 0, 0
    else:
        #pdb.set_trace()
        seq_array   = numpy.array(list(seq_match_span))
        tag_array   = numpy.array(list(tag_match_span))
        matches     = sum(seq_array == tag_array)
        error       = sum(seq_array != tag_array) + (len(tag) - \
            len(tag_match_span.replace('-','')))
        # I didn't like the way that the original method at
        # http://github.com/chapmanb/bcbb/tree/master treats gaps 
        return matches, error

def smithWaterman(seq, tags, allowed_errors):
    '''Smith-Waterman alignment method for aligning tags with their respective
    sequences.  Only called when regular expression matching patterns fail.
    Inspired by http://github.com/chapmanb/bcbb/tree/master'''
    high_score = {'tag':None, 'seq_match':None, 'mid_match':None, 'score':None, 
        'start':None, 'end':None, 'matches':None, 'errors':allowed_errors}
    for tag in tags:
        seq_match, tag_match, score, start, end = pairwise2.align.localms(seq, 
        tag, 5.0, -4.0, -9.0, -0.5, one_alignment_only=True)[0]
        seq_match_span  = seq_match[start:end]
        tag_match_span  = tag_match[start:end]
        match, errors   = matches(tag, seq_match_span, tag_match_span, allowed_errors)
        if match >= len(tag)-allowed_errors and match > high_score['matches'] \
            and errors <= high_score['errors']:
            high_score['tag'] = tag
            high_score['seq_match'] = seq_match
            high_score['tag_match'] = tag_match
            high_score['score'] = score
            high_score['start'] = start
            high_score['end'] = end
            high_score['matches'] = match
            high_score['seq_match_span'] = seq_match_span
            high_score['errors'] = errors
    if high_score['matches']:
        return high_score['tag'], high_score['matches'], \
        high_score['seq_match'], high_score['seq_match_span'], \
        high_score['start'], high_score['end']
    else:
        return None

def qualTrimming(sequence, min_score=10):
    '''Remove ambiguous bases from 5' and 3' sequence ends'''
    s = str(sequence.seq)
    sl = list(s)
    for q in enumerate(sequence.letter_annotations["phred_quality"]):
        if q[1] < min_score:
            sl[q[0]] = 'N'
    s = ''.join(sl)
    # find runs of ambiguous bases at 5' and 3' ends
    left_re, right_re = re.compile('^N+'),re.compile('N+$')
    left_trim, right_trim = re.search(left_re, s), re.search(right_re, s)
    if left_trim:
        left_trim = left_trim.end()
    if right_trim:
        right_trim = right_trim.end()
    return trim(sequence, left_trim, right_trim)

def midTrim(sequence, tags, max_gap_char=5, **kwargs):
    '''Remove the MID tag from the sequence read'''
    #if sequence.id == 'MID_No_Error_ATACGACGTA':
    #    pdb.set_trace()
    s = str(sequence.seq)
    mid = leftLinker(s, tags, max_gap_char, True, fuzzy=kwargs['fuzzy'])
    if mid:
        trimmed = trim(sequence, mid[3])
        tag, m_type, seq_match = mid[0], mid[1], mid[4]
        return tag, trimmed, seq_match, m_type
    else:
        return None

def SWMatchPos(seq_match_span, start, stop):
    # slice faster than ''.startswith()
    if seq_match_span[0] == '-':
        start = start + seq_match_span.count('-')
    else:
        stop = stop - seq_match_span.count('-')
    return start, stop

def leftLinker(s, tags, max_gap_char, gaps=False, **kwargs):
    '''Matching methods for left linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed'''
    for tag in tags:
        if gaps:
            r = re.compile(('^%s') % (tag))
        else:
            r = re.compile(('^[acgtnACGTN]{0,%s}%s') % (max_gap_char, tag))
        match = re.search(r, s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    #if s == 'ACCTCGTGCGGAATCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG':
    #    pdb.set_trace()
    if not match and kwargs['fuzzy']:
        match = smithWaterman(s, tags, 1)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return tag, m_type, start, stop, seq_match
    else:
        return None

def rightLinker(s, tags, max_gap_char, gaps=False, **kwargs):
    '''Mathing methods for right linker - regex first, followed by fuzzy (SW)
    alignment, if the option is passed'''
    #if s == 'GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG':
    #    pdb.set_trace()
    revtags = revCompTags(tags)
    for tag in revtags:
        if gaps:
            r = re.compile(('%s$') % (tag))
        else:
            r = re.compile(('%s[acgtnACGTN]{0,%s}$') % (tag, max_gap_char))
        match = re.search(r, s)
        if match:
            m_type = 'regex'
            start, stop = match.start(), match.end()
            # by default, this is true
            seq_match = tag
            break
    if not match and kwargs['fuzzy']:
        match = smithWaterman(s, revtags, 1)
        # we can trim w/o regex
        if match:
            m_type = 'fuzzy'
            tag = match[0]
            seq_match = match[3]
            start, stop = SWMatchPos(match[3],match[4], match[5])
    if match:
        return revComp(tag), m_type, start, stop, seq_match
    else:
        return None

def linkerTrim(sequence, tags, max_gap_char=22, **kwargs):
    '''Use regular expression and (optionally) fuzzy string matching
    to locate and trim linkers from sequences'''
    m_type  = False
    s       = str(sequence.seq)
    left    = leftLinker(s, tags, max_gap_char=22, fuzzy=kwargs['fuzzy'])
    right   = rightLinker(s, tags, max_gap_char=22, fuzzy=kwargs['fuzzy'])
    if left and right and left[0] == right[0]:
        # we can have lots of conditional matches here
        if left[2] <= max_gap_char and right[2] >= (len(s) - (len(right[0]) +\
        max_gap_char)):
            trimmed = trim(sequence, left[3], right[2])
            # left and right are identical so largely pass back the left
            # info... except for m_type which can be a combination
            tag, m_type, seq_match = left[0], left[1]+'-'+right[1]+'-both', \
            left[4]
        else:
            pass
    elif left and right and left[0] != right[0]:
        # flag
        if left[2] <= max_gap_char and right[2] >= (len(s) - (len(right[0]) +\
        max_gap_char)):
            trimmed = None
            tag, m_type, seq_match = None, 'tag-mismatch', None
    elif left:
        if left[2] <= max_gap_char:
            trimmed = trim(sequence, left[3])
            tag, m_type, seq_match = left[0], left[1]+'-left', left[4]
        else:
            # flag
            pass
    elif right:
        if right[2] >= (len(s) - (len(right[0]) + max_gap_char)):
            trimmed = trim(sequence, None, right[2])
            tag, m_type, seq_match = right[0], right[1]+'-right', right[4]
        else:
            # flag
            pass
    if m_type:
        try:
            return tag, trimmed, seq_match, tags[tag], m_type
        except:
            return tag, trimmed, seq_match, None, m_type
    else:
        return None

def reverse(items, null=False):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    if null:
        items += ((None, None),)
    for i in items:
        t = (i[1],i[0])
        l.append(t)
    return dict(l)

def createSeqTable(c):
    '''Create necessary tables in our database to hold the sequence and 
    tagging data'''
    # TODO:  move blob column to its own table, indexed by id
    # DONE:  move all tables to InnoDB??
    try:
        c.execute('''DROP TABLE sequence''')
    except:
        pass
    c.execute('''CREATE TABLE sequence (id INT UNSIGNED NOT NULL 
        AUTO_INCREMENT,name VARCHAR(100),mid VARCHAR(30),mid_seq VARCHAR(30),
        mid_match VARCHAR(30),mid_method VARCHAR(50),linker VARCHAR(50),
        linker_seq VARCHAR(50),linker_match VARCHAR(50),linker_method 
        VARCHAR(50),cluster VARCHAR(75),concat_seq VARCHAR(50), 
        concat_match varchar(50), concat_method VARCHAR(50),
        n_count SMALLINT UNSIGNED, untrimmed_len SMALLINT UNSIGNED, 
        seq_trimmed TEXT, trimmed_len SMALLINT UNSIGNED, record BLOB, PRIMARY
        KEY (id), INDEX sequence_cluster (cluster)) ENGINE=InnoDB''')

def createQualSeqTable(c):
    # TODO:  move blob column to its own table, indexed by id
    # DONE:  move all tables to InnoDB??
    try:
        c.execute('''DROP TABLE sequence''')
    except:
        pass
    c.execute('''CREATE TABLE sequence (id INT UNSIGNED NOT NULL 
        AUTO_INCREMENT,name VARCHAR(100), n_count SMALLINT UNSIGNED, 
        untrimmed_len MEDIUMINT UNSIGNED, seq_trimmed MEDIUMTEXT, trimmed_len 
        MEDIUMINT UNSIGNED, record MEDIUMBLOB, PRIMARY KEY (id)) ENGINE=InnoDB''')

def concatCheck(sequence, all_tags, all_tags_regex, reverse_linkers, **kwargs):
    '''Check screened sequence for the presence of concatemers by scanning 
    for all possible tags - after the 5' and 3' tags have been removed'''
    s = str(sequence.seq)
    m_type = None
    # do either/or to try and keep speed up, somewhat
    #if not kwargs['fuzzy']:
    #pdb.set_trace()
    for tag in all_tags_regex:
        match = re.search(tag, s)
        if match:
            tag = tag.pattern
            m_type = 'regex-concat'
            seq_match = tag
            break
    if not match and ['fuzzy']:
    #else:
        match = smithWaterman(s, all_tags, 1)
        # we can trim w/o regex
        if match:
            tag = match[0]
            m_type = 'fuzzy-concat'
            seq_match = match[3]
    if m_type:
        return tag, m_type, seq_match
    else:
        return None, None, None

def sequenceCount(input):
    '''Determine the number of sequence reads in the input'''
    handle = open(input, 'rU')
    lines = handle.read().count('>')
    handle.close()
    return lines
            

def qualOnlyWorker(sequence, qual, conf):
    # we need a separate connection for each mysql cursor or they are going
    # start going into locking hell and things will go poorly. Creating a new 
    # connection for each worker process is the easiest/laziest solution.
    # Connection pooling (DB-API) didn't work so hot, but probably because 
    # I'm slightly retarded.
    conn = MySQLdb.connect(user=conf.get('Database','USER'), 
        passwd=conf.get('Database','PASSWORD'), 
        db=conf.get('Database','DATABASE'))
    cur = conn.cursor()
    # convert low-scoring bases to 'N'
    untrimmed_len = len(sequence.seq)
    qual_trimmed = qualTrimming(sequence, qual)
    N_count = str(qual_trimmed.seq).count('N')
    sequence = qual_trimmed
    # pickle the sequence record, so we can store it as a BLOB in MySQL, we
    # can thus recurrect it as a sequence object when we need it next.
    sequence_pickle = cPickle.dumps(sequence,1)
    cur.execute('''INSERT INTO sequence (name, n_count, untrimmed_len, 
        seq_trimmed, trimmed_len, record) 
        VALUES (%s,%s,%s,%s,%s,%s)''', 
        (sequence.id, N_count, untrimmed_len, sequence.seq, len(sequence.seq), 
        sequence_pickle))
    #pdb.set_trace()
    cur.close()
    conn.commit()
    # keep our connection load low
    conn.close()
    return

class Record():
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, sequence):
        # super(Params, self).__init__()
        assert isinstance(sequence,SeqRecord), \
            'The Record class must be instantiated with a BioPython Seq object'
        self.unmod              = sequence # a biopython sequence object
        self.sequence           = None
        self.nCount             = None
        self.mid                = None
        self.mid_seq            = None
        self.reverse_mid        = None
        self.seq_match          = None
        self.m_type             = None
        self.l_seq              = None
        self.l_tag              = None
        self.l_seq_match        = None
        self.l_critter          = None
        self.l_m_type           = None
        self.reverse_linker     = None
        self.concat_tag         = None
        self.concat_type        = None
        self.concat_seq_type    = None
        self.concat_count       = None
        self.concat_seq_match   = None
    
    def __repr__(self):
        return '''<linkers.record for %s>''' % self.unmod.id

def linkerWorker(sequence, params):
    # we need a separate connection for each mysql cursor or they are going
    # start going into locking hell and things will go poorly. Creating a new 
    # connection for each worker process is the easiest/laziest solution.
    # Connection pooling (DB-API) didn't work so hot, but probably because 
    # I'm slightly retarded.
    conn = MySQLdb.connect(
        user=params.user,
        passwd=params.pwd,
        db=params.db
        )
    cur = conn.cursor()
    # for now, we'll keep this here
    seqRecord = Record(sequence)
    if params.qualTrim:
        seqRecord.sequence = qualTrimming(seqRecord.unmod, params.minQual)
    else:
        seqRecord.sequence = seqRecord.unmod
    seqRecord.nCount = str(seqRecord.sequence.seq).count('N')
    #pdb.set_trace()
    tags = params.tags
    if params.midTrim:
        # search on 5' (left) end for MID
        mid = midTrim(seqRecord.sequence, params.tags, params.midGap, fuzzy=params.fuzzy)
        if mid:
            # if MID, search for exact matches (for and revcomp) on Linker
            # provided no exact matches, use fuzzy matching (Smith-Waterman) +
            # error correction to find Linker
            seqRecord.mid           = mid[0]
            seqRecord.sequence      = mid[1]
            seqRecord.seq_match     = mid[2]
            seqRecord.m_type        = mid[3]
            seqRecord.reverse_mid   = params.reverse_mid[seqRecord.mid]
            tags                    = params.tags[seqRecord.mid]
    #pdb.set_trace()
    if params.linkerTrim:
        linker = linkerTrim(seqRecord.sequence, tags, params.linkerGap, fuzzy=params.fuzzy)
        if linker:
            if linker[0]:
                seqRecord.l_tag             = linker[0]
                seqRecord.sequence          = linker[1]
                seqRecord.l_seq_match       = linker[2]
                seqRecord.l_critter         = linker[3]
                seqRecord.l_m_type          = linker[4]
                seqRecord.reverse_linker    = params.reverse_linkers[seqRecord.l_tag]
            # deal with tag-mismatch
            if not linker[0] and linker[4]:
                seqRecord.l_m_type          = linker[4]
    # check for concatemers
    if params.concat:
        if l_trimmed and len(l_trimmed.seq) > 0:
            concat_tag, concat_type, concat_seq_match = concatCheck(l_trimmed, 
                all_tags, all_tags_regex, reverse_linkers, fuzzy=params.fuzzy)
        else:
            concat_tag, concat_type, concat_seq_match = None, None, None
    else:
        concat_tag, concat_type, concat_seq_match = None, None, None
    # pickle the sequence record, so we can store it as a BLOB in MySQL, we
    # can thus resurrect it as a sequence object when we need it next.
    sequence_pickle = cPickle.dumps(seqRecord.sequence,1)
    cur.execute('''INSERT INTO sequence (name, mid, mid_seq, mid_match, 
        mid_method, linker, linker_seq, linker_match, linker_method, cluster, 
        concat_seq, concat_match, concat_method, n_count, untrimmed_len, 
        seq_trimmed, trimmed_len, record) 
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''', 
        (seqRecord.sequence.id, seqRecord.reverse_mid, seqRecord.mid, \
        seqRecord.seq_match, seqRecord.m_type, seqRecord.reverse_linker, \
        seqRecord.l_tag, seqRecord.l_seq_match, seqRecord.l_m_type, \
        seqRecord.l_critter, seqRecord.concat_tag, \
        seqRecord.concat_seq_match, seqRecord.concat_type, seqRecord.nCount, \
        len(seqRecord.unmod.seq), seqRecord.sequence.seq, \
        len(seqRecord.sequence.seq), sequence_pickle))
    #pdb.set_trace()
    cur.close()
    conn.commit()
    # keep our connection load low
    conn.close()
    return



def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg

class Parameters():
    '''linkers.py run parameters'''
    def __init__(self, conf):
        self.conf            = conf
        self.db              = self.conf.get('Database','DATABASE')
        self.user            = self.conf.get('Database','USER')
        self.pwd             = self.conf.get('Database','PASSWORD')
        self.qualTrim        = self.conf.getboolean('Steps', 'Trim')
        self.minQual         = self.conf.getint('GeneralParameters', 'MinQualScore')
        self.midTrim         = self.conf.getboolean('Steps','MidTrim')
        self.midGap          = self.conf.getint('GeneralParameters','MidGap')
        self.linkerTrim      = self.conf.getboolean('Steps', 'LinkerTrim')
        self.linkerGap       = self.conf.getint('GeneralParameters','LinkerGap')
        self.concat          = self.conf.getboolean('GeneralParameters','CheckForConcatemers')
        self.fuzzy           = self.conf.getboolean('GeneralParameters','FuzzyMatching')
        self.mids            = None
        self.reverse_mid     = None
        self.linkers         = None
        self.reverse_linkers = None
        self.clust           = None
        self.tags            = None
        self.all_tags        = None
        self.all_tags_regex  = None
        self._setup()
    
    def __repr__(self):
        return '''<linkers.parameters run values>'''
    
    def _setup(self):
        if self.midTrim and self.linkerTrim:
            self._mid()
            self._linkers()
            self.clust       = self.conf.items('MidLinkerGroups')
            self._tagLibrary()

        elif self.midTrim and not self.LinkerTrim:
            self_mid()
            self.clust       = self.conf.items('MidGroup')
            self._tagLibrary()
            
        elif not self.midTrim and self.linkerTrim:
            self._linkers()
            self.clust       = self.conf.items('LinkerGroups')
            self._tagLibrary()
            
        # do we check for concatemers?
        if self.concat:
            self._allPossibleTags()
    
    def _linkers(self):
        self.linkers         = dict(self.conf.items('Linker'))
        self.reverse_linkers = reverse(self.conf.items('Linker'), True)
    
    def _mid(self):
        self.mids            = dict(self.conf.items('Mid'))
        self.reverse_mid     = reverse(self.conf.items('Mid'), True)
    
    def _tagLibrary(self):
        '''Create a tag-library from the mids and the linkers which allows us to 
        track which organisms go with which MID+linker combo'''
        self.tags = {}
        for c in self.clust:
            if self.mids and self.linkers:
                m,l = c[0].replace(' ','').split(',')
                org = c[1]
                if self.mids[m] not in self.tags.keys():
                    self.tags[self.mids[m]] = {self.linkers[l]:org}
                else:
                    self.tags[self.mids[m]][self.linkers[l]] = org
            
            elif not self.mids and self.linkers:
                l = c[0]
                org = c[1]
                self.tags[self.linkers[l]] = org
            
            elif self.mids and not self.linker:
                l = c[0]
                org = c[1]
                self.tags[self.mids[l]] = org
    
    def _allPossibleTags(self):
        '''Create regular expressions for the forward and reverse complements
        of all of the tags sequences used in a run'''
        # at = all tags; rat = reverse complement all tags
        self.all_tags = []
        self.all_tags_regex = []
        for c in self.clust:
            if self.mids and self.linkers:
                m,l = c[0].replace(' ','').split(',')
            elif not self.mids and self.linkers:
                l = c[0]
            elif self.mids and not self.linkers:
                pass
            all_tags.append(self.linkers[l])
            all_tags_regex.append(re.compile('%s' % self.linkers[l]))
            all_tags.append(revComp(self.linkers[l]))
            all_tags_regex.append(re.compile('%s' % revComp(self.linkers[l])))

def main():
    '''Main loop'''
    start_time      = time.time()
    options, arg    = interface()
    motd()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    conf            = ConfigParser.ConfigParser()
    conf.read(options.conf)
    # build our configuration
    params = Parameters(conf)
    conn = MySQLdb.connect(
        user=params.user,
        passwd=params.pwd,
        db=params.db
        )
    cur = conn.cursor()
    # crank out a new table for the data
    createSeqTable(cur)
    conn.commit()
    seqcount = sequenceCount(conf.get('Input','sequence'))
    sequence = QualityIO.PairedFastaQualIterator(
    open(conf.get('Input','sequence'), "rU"), 
    open(conf.get('Input','qual'), "rU"))
    #pdb.set_trace()
    if conf.getboolean('Multiprocessing', 'MULTIPROCESSING'):
        # get num processors
        n_procs = conf.get('Multiprocessing','processors')
        if n_procs == 'Auto':
            # we'll use x-1 cores (where x = avail. cores)
            n_procs = multiprocessing.cpu_count() - 1
        else:
            n_procs = int(n_procs)
        print 'Multiprocessing.  Number of processors = ', n_procs
        # to test with fewer sequences
        #count = 0
        try:
            threads = []
            pb = progress.bar(0,seqcount,60)
            pb_inc = 0
            while sequence:
                if len(threads) < n_procs:
                    p = multiprocessing.Process(
                            target=linkerWorker, 
                            args=(
                                sequence.next(),
                                params,
                                )
                            )
                    p.start()
                    threads.append(p)
                    if (pb_inc+1)%1000 == 0:
                        pb.__call__(pb_inc+1)
                    elif pb_inc + 1 == seqcount:
                        pb.__call__(pb_inc+1)
                    pb_inc += 1
                else:
                    for t in threads:
                        if not t.is_alive():
                            threads.remove(t)
        except StopIteration:
            pass
    else:
        print 'Not using multiprocessing'
        count = 0
        try:
            pb = progress.bar(0,seqcount,60)
            pb_inc = 0
            #while count < 1000:
            while sequence:
                #count +=1
                linkerWorker(sequence.next(), params)
                if (pb_inc+1)%1000 == 0:
                    pb.__call__(pb_inc+1)
                elif pb_inc + 1 == seqcount:
                    pb.__call__(pb_inc+1)
                pb_inc += 1
        except StopIteration:
            pass
    print '\n'
    cur.close()
    conn.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print '\nTime for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()
