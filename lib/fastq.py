#!/usr/bin/env python
# encoding: utf-8

"""
fastq.py

These files (fasta.py, fastq.py, transform.py, sequence.py) are originally part of the galaxy-dist 
package (http://bitbucket.org/galaxy/galaxy-dist/src/), and these particular classes, methods, and 
functions were created by Dan Blankenberg.  

See Blankenberg et al. doi:  10.1093/bioinformatics/btq281

I have modified the original source, changing the way that quality scores are stored and used
(all quality scores as numpy arrays in sanger phred values with methods to convert and display as
string), and addind additional method to several of the classes (e.g. trimming within 
sequence.SequenceRead - subclassed by both fasta and fastq).  I also altered methods and class 
method to use numpy methods within functions, where possible.  Finally, i slightly changed how 
how fasta and fastq files are read.

The original the code was Copyright (c) 2005 Pennsylvania State University.  See my 
LICENSE file for the license and Copyright of the current code (also BSD).

This file incorporates code covered by the following terms:

        Copyright (c) 2005 Pennsylvania State University

        Permission is hereby granted, free of charge, to any person obtaining
        a copy of this software and associated documentation files (the
        "Software"), to deal in the Software without restriction, including
        without limitation the rights to use, copy, modify, merge, publish,
        distribute, sublicense, and/or sell copies of the Software, and to
        permit persons to whom the Software is furnished to do so, subject to
        the following conditions:

        The above copyright notice and this permission notice shall be
        included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
        EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
        IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
        CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
        SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

import math
import copy
import numpy
import string
import transform
from sequence import SequencingRead

import pdb

class fastqSequencingRead(SequencingRead):
    format = 'sanger'
    ascii_min = 33
    ascii_max = 126
    quality_min = 0
    quality_max = 93
    score_system = 'phred'
    sequence_space = 'base'
    
    @classmethod
    def get_class_by_format(cls, format):
        """return class of a format"""
        assert format in FASTQ_FORMATS.keys(), 'Unknown format type specified: %s' % format
        return FASTQ_FORMATS[ format ]
    
    @classmethod
    def convert_score_phred_to_solexa(cls, quality_array):
        """convert phred quality scores to solexa quality scores and return array of values"""
        quality_array[numpy.where(quality_array <= 0)[0]] = -5
        return numpy.around(10.0 * numpy.log10(pow(10.0, quality_array/10.0) - 1.0)).astype('uint8')
    
    @classmethod
    def convert_score_solexa_to_phred(cls, quality_array):
        """convert solexa quality scores to phred quality scores and return array of values"""
        return numpy.around(10.0 * numpy.log10(pow(10.0, quality_array/10.0) + 1.0)).astype('uint8')
    
    @classmethod
    def restrict_scores_to_valid_range(cls, quality_array):
        """limit scores in the quality array to the class maximum and return ceil(array)"""
        return quality_array.clip(max=cls.quality_max)

    def set_quality(self, quality_string):
        """parse ASCII and QUAL strings into and unsigned integer array and set the quality 
        attribute.  typically overrides sequence.SequencingRead.set_quality()"""
        quality_string = quality_string.strip()
        if ' ' in quality_string or '\t' in quality_string:
            self.quality = self._get_qual_array_from_decimal(quality_string)
        elif quality_string:
            temp_array = numpy.fromstring(quality_string, dtype = 'uint8')
            assert self.is_valid_format(temp_array), "Ascii quality string is not valid, given format."
            self.quality = self._get_qual_array_from_ascii(temp_array)
            self.assert_sequence_quality_lengths()
        else:
            self.quality = None
    
    def _get_qual_array_from_ascii(self, quality_array):
        """PRIVATE helper to parse an ASCII quality string and return an unsigned integer array of quality values"""
        return (quality_array - self.ascii_min + self.quality_min).astype('uint8')
    
    def get_quality_string(self, qtype = 'decimal'):
        """format and return the quality array as string in decimal (default) or ascii format"""
        if self.quality.any() and qtype == 'decimal':
            return ' '.join(self.quality.astype('|S2'))
        elif self.quality.any() and qtype == 'ascii':
            # numpy can return null bytes when there on non-letter characters in an array (e.g. 
            # "@".  strip those, return the string.)
            return (self.quality + self.ascii_min - (self.quality_min)).astype('uint8').tostring()
        else:
            return None

    def convert_read_to_format(self, format):
        """convert reads between formats and return a new fastq object (of `format` class)"""
        assert format in FASTQ_FORMATS, 'Unknown format type specified: %s' % format
        new_class = FASTQ_FORMATS[ format ]
        new_read = new_class()
        new_read.identifier = self.identifier
        # we've removed color space
        new_read.sequence = self.sequence
        new_read.description = self.description
        if self.score_system == 'phred' and new_read.score_system == 'solexa':
            temp = self.convert_score_phred_to_solexa(self.quality)
        elif self.score_system == 'solexa' and new_read.score_system == 'phred':
            temp = self.convert_score_solexa_to_phred(self.quality)
        else:
            temp = self.quality
        new_read.quality = new_class.restrict_scores_to_valid_range(temp)
        return new_read
    
    def get_sequence(self):
        """return sequence"""
        return self.sequence
    
    #def slice(self, left_column_offset, right_column_offset):
    #    """slice sequence (and quality) from left, right, or left and right and returna new fastqSequencingRead object"""
    #    new_read = fastqSequencingRead.get_class_by_format(self.format)()
    #    new_read.identifier = self.identifier
    #    new_read.sequence = self.get_sequence()[left_column_offset:right_column_offset]
    #    new_read.description = self.description
    #    new_read.quality = self.quality[left_column_offset:right_column_offset]
    #    return new_read
    
    def is_valid_format(self, quality_string):
        """ensure that quality values are within min-max range and return boolean result"""
        if ((self.ascii_min <= quality_string) & (quality_string <= self.ascii_max)).all():
            return True
        else:
            return False
        #if self.is_ascii_encoded():
        #    for val in self.get_ascii_quality_scores():
        #        val = ord( val )
        #        if val < self.ascii_min or val > self.ascii_max:
        #            return False
        #else:
        #    for val in self.get_decimal_quality_scores():
        #        if val < self.quality_min or val > self.quality_max:
        #            return False
        #if not self.is_valid_sequence():
        #    return False
        #return True
    
    def is_valid_sequence( self ):
        for base in self.get_sequence():
            if base not in self.valid_sequence_list:
                return False
        return True
    
    def insufficient_quality_length( self ):
        return len(self.quality) < len(self.sequence)
    
    def assert_sequence_quality_lengths(self):
        q,s = len(self.quality), len(self.sequence)
        assert s == q, "Invalid FASTQ file: quality score length ({0}) does not match sequence length ({1})".format(q,s)
    
    def reverse( self, clone = True ):
        # TODO:  remove this method?
        #need to override how decimal quality scores are reversed
        if clone:
            rval = self.clone()
        else:
            rval = self
        rval.sequence = transform.reverse( self.sequence )
        if rval.is_ascii_encoded():
            rval.quality = rval.quality[::-1]
        else:
            rval.quality = reversed( rval.get_decimal_quality_scores() )
            rval.quality = "%s " % " ".join( map( str, rval.quality ) )
        return rval

class fastqSangerRead( fastqSequencingRead ):
    format = 'sanger'
    ascii_min = 33
    ascii_max = 126
    quality_min = 0
    quality_max = 93
    score_system = 'phred'
    sequence_space = 'base'

class fastqIlluminaRead( fastqSequencingRead ):
    format = 'illumina'
    ascii_min = 64
    ascii_max = 126
    quality_min = 0
    quality_max = 62
    score_system = 'phred'
    sequence_space = 'base'

class fastqSolexaRead( fastqSequencingRead ):
    format = 'solexa'
    ascii_min = 59
    ascii_max = 126
    quality_min = -5
    quality_max = 62
    score_system = 'solexa'
    sequence_space = 'base'

FASTQ_FORMATS = {}
for format in [ fastqIlluminaRead, fastqSolexaRead, fastqSangerRead ]:
    FASTQ_FORMATS[ format.format ] = format


class fastqAggregator( object ):
    VALID_FORMATS = FASTQ_FORMATS.keys()
    def __init__( self,  ):
        self.ascii_values_used = [] #quick lookup of all ascii chars used
        self.seq_lens = {} #counts of seqs by read len
        self.nuc_index_quality = [] #counts of scores by read column
        self.nuc_index_base = [] #counts of bases by read column
    def consume_read( self, fastq_read ):
        #ascii values used
        for val in fastq_read.get_ascii_quality_scores():
            if val not in self.ascii_values_used:
                self.ascii_values_used.append( val )
        #lengths
        seq_len = len( fastq_read )
        self.seq_lens[ seq_len ] = self.seq_lens.get( seq_len, 0 ) + 1
        #decimal qualities by column
        for i, val in enumerate( fastq_read.get_decimal_quality_scores() ):
            if i == len( self.nuc_index_quality ):
                self.nuc_index_quality.append( {} )
            self.nuc_index_quality[ i ][ val ] = self.nuc_index_quality[ i ].get( val, 0 ) + 1
        #bases by column
        for i, nuc in enumerate( fastq_read.get_sequence() ):
            if i == len( self.nuc_index_base ):
                self.nuc_index_base.append( {} )
            nuc = nuc.upper()
            self.nuc_index_base[ i ][ nuc ] = self.nuc_index_base[ i ].get( nuc, 0 ) + 1
    def get_valid_formats( self, check_list = None ):
        if not check_list:
            check_list = self.VALID_FORMATS
        rval = []
        sequence = []
        for nuc_dict in self.nuc_index_base:
            for nuc in nuc_dict.keys():
                if nuc not in sequence:
                    sequence.append( nuc )
        sequence = "".join( sequence )
        quality = "".join( self.ascii_values_used )
        for fastq_format in check_list:
            fastq_read = fastqSequencingRead.get_class_by_format( fastq_format )()
            fastq_read.quality = quality
            fastq_read.sequence = sequence
            if fastq_read.is_valid_format():
                rval.append( fastq_format )
        return rval
    def get_ascii_range( self ):
        return ( min( self.ascii_values_used ), max( self.ascii_values_used ) )
    def get_decimal_range( self ):
        decimal_values_used = []
        for scores in self.nuc_index_quality:
            decimal_values_used.extend( scores.keys() )
        return ( min( decimal_values_used ), max( decimal_values_used ) )
    def get_length_counts( self ):
        return self.seq_lens
    def get_max_read_length( self ):
        return len( self.nuc_index_quality )
    def get_read_count_for_column( self, column ):
        if column >= len( self.nuc_index_quality ):
            return 0
        return sum( self.nuc_index_quality[ column ].values() )
    def get_read_count( self ):
        return self.get_read_count_for_column( 0 )
    def get_base_counts_for_column( self, column ):
        return self.nuc_index_base[ column ]
    def get_score_list_for_column( self, column ):
        return self.nuc_index_quality[ column ].keys()
    def get_score_min_for_column( self, column ):
        return min( self.nuc_index_quality[ column ].keys() )
    def get_score_max_for_column( self, column ):
        return max( self.nuc_index_quality[ column ].keys() )
    def get_score_sum_for_column( self, column ):
        return sum( score * count for score, count in self.nuc_index_quality[ column ].iteritems() )
    def get_score_at_position_for_column( self, column, position ):
        score_value_dict = self.nuc_index_quality[ column ]
        scores = sorted( score_value_dict.keys() )
        for score in scores:
            if score_value_dict[ score ] <= position:
                position -= score_value_dict[ score ]
            else:
                return score
    def get_summary_statistics_for_column( self, i ):
        def _get_med_pos( size ):
            halfed = int( size / 2 )
            if size % 2 == 1:
                return [ halfed ]
            return[ halfed - 1, halfed ]
        read_count = self.get_read_count_for_column( i )
        
        min_score = self.get_score_min_for_column( i )
        max_score = self.get_score_max_for_column( i )
        sum_score = self.get_score_sum_for_column( i )
        mean_score = float( sum_score ) / float( read_count )
        #get positions
        med_pos = _get_med_pos( read_count )
        if 0 in med_pos:
            q1_pos = [ 0 ]
            q3_pos = [ read_count - 1 ]
        else:
            q1_pos = _get_med_pos( min( med_pos ) )
            q3_pos = []
            for pos in q1_pos:
                q3_pos.append( max( med_pos ) + 1 + pos )
        #get scores at position
        med_score = float( sum( [ self.get_score_at_position_for_column( i, pos ) for pos in med_pos ] ) ) / float( len( med_pos ) )
        q1 = float( sum( [ self.get_score_at_position_for_column( i, pos ) for pos in q1_pos ] ) ) / float( len( q1_pos ) )
        q3 = float( sum( [ self.get_score_at_position_for_column( i, pos ) for pos in q3_pos ] ) ) / float( len( q3_pos ) )
        #determine iqr and step
        iqr = q3 - q1
        step = 1.5 * iqr
        
        #Determine whiskers and outliers
        outliers = []
        score_list = sorted( self.get_score_list_for_column( i ) )
        left_whisker = q1 - step
        for score in score_list:
            if left_whisker <= score:
                left_whisker = score
                break
            else:
                outliers.append( score )
        
        right_whisker = q3 + step
        score_list.reverse()
        for score in score_list:
            if right_whisker >= score:
                right_whisker = score
                break
            else:
                outliers.append( score )
        
        column_stats = { 'read_count': read_count, 
                         'min_score': min_score, 
                         'max_score': max_score, 
                         'sum_score': sum_score, 
                         'mean_score': mean_score, 
                         'q1': q1, 
                         'med_score': med_score, 
                         'q3': q3, 
                         'iqr': iqr, 
                         'left_whisker': left_whisker, 
                         'right_whisker': right_whisker,
                         'outliers': outliers }
        return column_stats

class fastqReader( object ):
    
    def __init__( self, fh, format = 'sanger' ):
        self.file = open(fh)
        self.format = format
        
    def close( self ):
        return self.file.close()
        
    def next(self):
        while True:
            fastq_header = self.file.readline()
            if not fastq_header:
                raise StopIteration
            fastq_header = fastq_header.rstrip( '\n\r' )
            #remove empty lines, apparently extra new lines at end of file is common?
            if fastq_header:
                break
                
        assert fastq_header.startswith( '@' ), 'Invalid fastq header: %s' % fastq_header
        rval = fastqSequencingRead.get_class_by_format(self.format)()        
        rval.identifier = fastq_header
        quality_string = ''
        while True:
            line = self.file.readline()
            if not line:
                raise Exception( 'Invalid FASTQ file: could not parse second instance of sequence identifier.' )
            line = line.rstrip( '\n\r' )
            if line.startswith( '+' ) and ( len( line ) == 1 or line[1:].startswith( fastq_header[1:] ) ):
                rval.description = line
                break
            rval.append_sequence( line )
        while len(quality_string) < len(rval.sequence):
            line = self.file.readline()
            if not line:
                break
            quality_string = ''.join([quality_string, line.strip()])
        rval.set_quality(quality_string)
        # ensure quality == sequence length
        rval.assert_sequence_quality_lengths()
        return rval
    
    def __iter__( self ):
        while True:
            yield self.next()

class fastqWriter( object ):
    
    def __init__( self, fh, format = None, force_quality_encoding = None ):
        self.file = fh
        self.format = format
        self.force_quality_encoding = force_quality_encoding
    
    def write( self, fastq_read ):
        if self.format:
            fastq_read = fastq_read.convert_read_to_format( self.format, force_quality_encoding = self.force_quality_encoding )
        self.file.write(str(fastq_read))
    
    def close( self ):
        return self.file.close()

class fastqJoiner( object ):
    
    def __init__( self, format, force_quality_encoding = None ):
        self.format = format
        self.force_quality_encoding = force_quality_encoding
    
    def join( self, read1, read2 ):
        if read1.identifier.endswith( '/2' ) and read2.identifier.endswith( '/1' ):
            #swap 1 and 2
            tmp = read1
            read1 = read2
            read2 = tmp
            del tmp
        if read1.identifier.endswith( '/1' ) and read2.identifier.endswith( '/2' ):
            identifier = read1.identifier[:-2]
        else:
            identifier = read1.identifier
        
        #use force quality encoding, if not present force to encoding of first read
        force_quality_encoding = self.force_quality_encoding
        if not force_quality_encoding:
            if read1.is_ascii_encoded():
                force_quality_encoding = 'ascii'
            else:
                force_quality_encoding = 'decimal'
        
        new_read1 = read1.convert_read_to_format( self.format, force_quality_encoding = force_quality_encoding )
        new_read2 = read2.convert_read_to_format( self.format, force_quality_encoding = force_quality_encoding )
        rval = FASTQ_FORMATS[ self.format ]()
        rval.identifier = identifier
        if len( read1.description ) > 1:
            rval.description = "+%s" % ( identifier[1:] )
        else:
            rval.description = '+'
        if rval.sequence_space == 'color':
            #need to handle color space joining differently
            #convert to nuc space, join, then convert back
            rval.sequence = rval.convert_base_to_color_space( new_read1.convert_color_to_base_space( new_read1.sequence ) + new_read2.convert_color_to_base_space( new_read2.sequence ) )
        else:
            rval.sequence = new_read1.sequence + new_read2.sequence
        if force_quality_encoding == 'ascii':
            rval.quality = new_read1.quality + new_read2.quality
        else:
            rval.quality = "%s %s" % ( new_read1.quality.strip(), new_read2.quality.strip() )
        return rval
    
    def get_paired_identifier( self, fastq_read ):
        identifier = fastq_read.identifier
        if identifier[-2] == '/':
            if identifier[-1] == "1":
                identifier = "%s2" % identifier[:-1]
            elif identifier[-1] == "2":
                identifier = "%s1" % identifier[:-1]
        return identifier

class fastqSplitter( object ):   
    def split( self, fastq_read ):
        length = len( fastq_read )
        #Only reads of even lengths can be split
        if length % 2 != 0:
            return None, None
        half = int( length / 2 )
        read1 = fastq_read.slice( 0, half )
        read1.identifier += "/1"
        if len( read1.description ) > 1:
            read1.description += "/1"
        read2 = fastq_read.slice( half, None )
        read2.identifier += "/2"
        if len( read2.description ) > 1:
            read2.description += "/2"
        return read1, read2
