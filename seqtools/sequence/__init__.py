"""

Dependencies:  numpy

"""


__author__ = 'Brant Faircloth'
__copyright__ = 'Copyright (c) 2010-2011, Brant C. Faircloth'
__credits__ = ['Brant Faircloth']
__license__ = 'BSD'
__version__ = '0.1'
__maintainer__ = 'Brant Faircloth'
__email__ = 'brant dot faircloth at gmail dot com'
__status__ = 'testing'

import sys

if sys.version_info < (2, 4):
    raise "must use python 2.5 or greater" # pragma: no cover

try: # pragma: no cover
    import numpy 
except ImportError: # pragma: no cover
    raise ImportError, 'NumPy does not seem to be installed. Please see the user guide.' # pragma: no cover

try: # pragma: no cover
    import multiprocessing 
except ImportError: # pragma: no cover
    raise ImportError, ''.join(['multiprocessing does not appear to be installed - either install ',
        'Python 2.6.x or greater, or try the multiprocessing backport (untested!):  ',
        'http://code.google.com/p/python-multiprocessing/'])
        
import fasta
import fastq
import fops
import sequence
import transform