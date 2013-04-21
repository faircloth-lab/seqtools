import os
import glob
import unittest

import pdb

def test():
    dirname = os.path.dirname(os.path.abspath(__file__))
    test_file_strings = glob.glob(os.path.join(dirname, 'test_*.py'))
    module_strings = [os.path.splitext(os.path.basename(name))[0]
            for name in test_file_strings]
    suites = [unittest.defaultTestLoader.loadTestsFromName(modstr) for modstr
            in module_strings]
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    test_suite = test()
    text_runner = unittest.TextTestRunner().run(test_suite)
