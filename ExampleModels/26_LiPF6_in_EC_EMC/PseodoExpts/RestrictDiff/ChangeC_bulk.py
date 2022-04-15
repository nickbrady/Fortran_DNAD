#
import sys
import numpy as np                          # Useful for numerical calculations on an array
import os                                   # Used for changing directories and open files with folders
from subprocess import call                 # Allows Python to call shell script
import shutil  								# allows python to execute terminal commands (such as rm mv cp ...)
import re 									# module used in sed (function) to find and replace text
import glob 								# module used to find files in a directory (used to remove temporary files)
from operator import itemgetter   			# Allows us to find groups of consecutive numbers
from tqdm import tqdm 						# module to produce and view a progress bar; valuable when running many simulations
import tempfile as tp
# import time

def sed(pattern, replace, source, dest=None, count=0):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count   (int):   number of occurrences to replace
        dest    (str):    destination filename, if not given, source will be over written.
    """

    fin = open(source, 'r')
    num_replaced = count

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = tp.mkstemp()
        fout = open(name, 'w')

    for line in fin:
        out = re.sub(pattern, replace, line)
        fout.write(out)

        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()

    if not dest:
        shutil.move(name, source)




cwd = os.getcwd()

_source = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_LiPF6_in_EC_EMC/With_Current/LiPF6_in_EC_EMC.f95'

destination = 'LiPF6_in_EC_EMC_val.f95'


run_fortran_path = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/RunFortran.py'

concentrations = [0.5, 1.0, 2.5]
for conc in concentrations:
    os.chdir(cwd)
    subfolder = '{}'.format(conc)

    try:
        os.makedirs(subfolder)
    except:
        pass

    os.chdir(subfolder)

    _destination = os.path.join(cwd, subfolder, destination)

    shutil.copyfile(_source, _destination)

    # need lower the applied current for lower concentrations
    parameter_strings = ['i_app_cm2 = 1e-4\\*8']
    sed(parameter_strings[0], 'i_app_cm2 = 1e-4*3', _destination)

    parameter_strings = ['c_LiPF6_bulk = 1.0e-3']
    sed(parameter_strings[0], 'c_LiPF6_bulk = {}e-3'.format(conc), _destination)


    p = call(['python3', run_fortran_path, destination])
