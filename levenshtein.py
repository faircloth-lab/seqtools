import pdb
import operator
import optparse
import os
import ConfigParser

def hamming_distance(s1, s2):
    '''Find the Hamming distance btw. 2 strings.
    
    Substitutions only.
    
    From http://en.wikipedia.org/wiki/Hamming_distance
    
    '''
    assert len(s1) == len(s2)
    return sum([ch1 != ch2 for ch1, ch2 in zip(s1, s2)])

def levenshtein(a,b):
    '''Calculates the levenshtein distance between a and b.
    
    Insertions, deletions, substitutions
    
    From here:  http://hetland.org/coding/python/levenshtein.py
    
    '''
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
    return current[n]

def getDistance(linkers, *args):
    linker_dist = []
    for l1 in xrange(len(linkers)):
        s1 = linkers[l1]
        for s2 in linkers[l1+1:]:
            edit_distance = levenshtein(s1[1],s2[1])
            linker_dist.append((s1[0], s2[0], edit_distance))
    #pdb.set_trace()
    link_list = [i[0] for i in linkers]
    min_dist = min([i[2] for i in linker_dist])
    if args[0] == 'midrl5':
        linker_sorted = sorted(linker_dist, key=operator.itemgetter(2))
        for item in linker_sorted:
            print item[0],item[1],item[2]
    return link_list, min_dist

def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', '-c', dest = 'conf', action='store', \
        type='string', default = None, \
        help='The path to the configuration file.', \
        metavar='FILE')
    p.add_option('--section', '-s', dest = 'section', action = 'store',\
        type='string', default = None, \
        help='The section of the config file to evaluate', metavar='FILE')
    (options,arg) = p.parse_args()
    if not options.conf:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.conf):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg 

def test():
    assert(levenshtein('shit','shat') == 1)
    assert(levenshtein('cool','bowl') == 2)
    assert(levenshtein('kitten','sitting') == 3)
    assert(levenshtein('bonjour','bougeoir') == 4)

def main():
    options, arg = interface()
    conf = ConfigParser.ConfigParser()
    conf.read(options.conf)
    groups = {}
    if options.section == 'Clusters':
        clust = conf.items('Clusters')
        links = dict(conf.items('Linker'))
        for c in clust:
            group = c[0].split(',')[0]
            linker = c[0].split(',')[1].strip()
            if group in groups.keys():
                groups[group] = groups[group] + ((linker,links[linker]),)
            else:
                groups[group] = ((linker,links[linker]),) 
    elif options.section == 'LinkerGroups':
        #pdb.set_trace()
        clust = conf.items('LinkerGroups')
        links = dict(conf.items('Linker'))
        g = ()
        for c in clust:
            g += ((c[0], links[c[0]]),)
        groups[1] = g
    for g in groups:
        #pdb.set_trace()
        ed = getDistance(groups[g], g)
        ed[0].sort()
        print '%s %s\n\tlinkers = %s\n\tedit distance = %s' % (options.section, g, str(ed[0]), ed[1])

    
if __name__ == '__main__':
    main()
    

'''linkers = (("SIMPLEX3", "CGTGCTGCGGAACT"),
("SIMPLEX6", "GCACGAGCGGAACT"),
("SIMPLEX11", "CGAGCGAGCGAAGT"),
("SIMPLEX12", "GCTGGCGTCGAAGT"),
("SIMPLEX13", "CCAGCACCGGAACA"),
("SIMPLEX14", "CCTGGGCACGAAGA"),
("SIMPLEX16", "GCAGCGTCGGAAAG"),
("SIMPLEX17", "GCTCCTGGCGAATC"))
original_run_linkers= (('SimpleX1','ACGTCGTGCGGAATC'), ('SimpleX2','AGCTGCTGGCGAATC'), ('SimpleX3','ACGTGCTGCGGAACT'), ('SimpleX4','AGCAGCAGCGGAATC'), ('SimpleX5','ACGAGCAGCGGAACT'), ('SimpleX6','AGCACGAGCGGAACT'), ('SuperSNX','GTTTAAGGCCTAGCTAGCAGAATC'))
getDistance(original_run_linkers)
MID = (('MID1', 'ACGAGTGCGT'), ('MID2', 'ACGCTCGACA'), ('MID3', 'AGACGCACTC'), ('MID4', 'AGCACTGTAG'), ('MID5', 'ATCAGACACG'), ('MID6', 'ATATCGCGAG'), ('MID7', 'CGTGTCTCTA'), ('MID8', 'CTCGCGTGTC'), ('MID10', 'TCTCTATGCG'), ('MID11', 'TGATACGTCT'))
getDistance(MID)'''

#Full Set
    
'''linkers = (("SIMPLEX3", "CGTGCTGCGGAACT"),
("SIMPLEX4", "GCAGCAGCGGAATC"),
("SIMPLEX6", "GCACGAGCGGAACT"),
("SIMPLEX7", "GGTCGAGCGGAATG"),
("SIMPLEX8", "GGTGCAGGCGAATG"),
("SIMPLEX11", "CGAGCGAGCGAAGT"),
("SIMPLEX12", "GCTGGCGTCGAAGT"),
("SIMPLEX13", "CCAGCACCGGAACA"),
("SIMPLEX14", "CCTGGGCACGAAGA"),
("SIMPLEX15", "CGTCGTGCGGAAAC"),
("SIMPLEX16", "GCAGCGTCGGAAAG"),
("SIMPLEX17", "GCTCCTGGCGAATC"))'''


#OLD Linkers

'''original_run_linkers= (('SimpleX1','ACGTCGTGCGGAATC'), ('SimpleX2','AGCTGCTGGCGAATC'), ('SimpleX3','ACGTGCTGCGGAACT'), ('SimpleX4','AGCAGCAGCGGAATC'), ('SimpleX5','ACGAGCAGCGGAACT'), ('SimpleX6','AGCACGAGCGGAACT'), ('SuperSNX','GTTTAAGGCCTAGCTAGCAGAATC'))'''
'''linkers = (('SimpleXL1', 'CGTCGTGCGGAATC'), ('SimpleXL2', 'GCTGCTGGCGAATC'), ('SimpleXL3', 'CGTGCTGCGGAACT'), ('SimpleXL4', 'GCAGCAGCGGAATC'), ('SimpleXL5', 'CGAGCAGCGGAACT'), ('SimpleXL6', 'GCACGAGCGGAACT'), ('SimpleXL7', 'GGTCGAGCGGAATG'), ('SimpleXL8', 'GGTGCAGGCGAATG'), ('SimpleXL9', 'CGTGCAGCGGAAGT'), ('SimpleXL10', 'GCAGCGTCGGAATG'), ('SimpleXL11', 'CGAGCGAGCGAAGT'), ('SimpleXL12', 'GCTGGCGTCGAAGT'), ('SimpleXL13', 'CCAGCAGCGGAACA'), ('SimpleXL14', 'CCTGCACGCGAAGA'), ('SimpleXL15', 'CGTCGTGCGGAAaC'), ('SimpleXL16', 'CGAGCGAGCGAAGT'))
'''

'''linkers = (('SimpleXL3', 'CGTGCTGCGGAACT'), ('SimpleXL4', 'GCAGCAGCGGAATC'), ('SimpleXL6', 'GCACGAGCGGAACT'), ('SimpleXL7', 'GGTCGAGCGGAATG'), ('SimpleXL8', 'GGTGCAGGCGAATG'), ('SimpleXL11', 'CGAGCGAGCGAAGT'), ('SimpleXL12', 'GCTGGCGTCGAAGT'), ('SimpleXL14', 'CCTGCACGCGAAGA'), ('SimpleXL15', 'CGTCGTGCGGAAAC'), ('SimpleXL16', 'GCAGCGTCGGAAaG'), ('SimpleXL17', 'GCTcCTGGCGAATC'), ('SimpleXL13b', 'CcAGCAcCGGAACa'))'''
