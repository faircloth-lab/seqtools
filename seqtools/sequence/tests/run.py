import glob
import unittest

def test():
    test_file_strings = glob.glob('test_*.py')
    module_strings = [str[0:len(str)-3] for str in test_file_strings]
    suites = [unittest.defaultTestLoader.loadTestsFromName(str) for str
            in module_strings]
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    test_suite = test()
    text_runner = unittest.TextTestRunner().run(test_suite)
