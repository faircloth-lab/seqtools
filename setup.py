import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup
from setuptools import find_packages

if __name__ == '__main__':
    setup(
        name='seqtools',
        version="0.5",
        description="Tools for sequence and alignment manipulation",
        author="Brant Faircloth",
        author_email="brant.faircloth+seqtools@gmail.com ",
        url="http://github.com/faircloth-lab/seqtools/",
        license="http://www.opensource.org/licenses/BSD-3-Clause",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        requires=['NumPy (>=1.3)',],
        long_description=open('README.txt').read(),
        package_data = {
                # If any package contains *.txt or *.rst files, include them:
                '': ['*.txt', '*.rst'],
                'sequence': ['tests/test-data/*'],
            },
        packages=[
                    'align',
                    'fs',
                    'sequence',
                    'sequence.tests'
                ],
        )
