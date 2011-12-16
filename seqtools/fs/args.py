
import os
import argparse


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def is_dir(dirname):
    """Using argparse, test to ensure something is a dir"""
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_sheet(filename):
    """Using argparse, test to ensure something is xlsx or csv"""
    if not os.path.basename(filename).endswith('csv') or not \
        os.path.basename(filename).endswith('xlsx'):
            msg = "{0} is not a sheet".format(filename)
            raise argparse.ArgumentTypeError(msg)
    else:
        return filename

