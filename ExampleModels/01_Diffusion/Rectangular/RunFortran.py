#
import numpy as np                          # Useful for numerical calculations on an array
import os                                   # Used for changing directories and open files with folders
from subprocess import call                 # Allows Python to call shell script
import shutil  								# allows python to execute terminal commands (such as rm mv cp ...)
import re 									# module used in sed (function) to find and replace text
import glob 								# module used to find files in a directory (used to remove temporary files)
from operator import itemgetter   			# Allows us to find groups of consecutive numbers
# from tqdm import tqdm 						# module to produce and view a progress bar; valuable when running many simulations
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

def run_fortran(_run_file_=None):
    # p = call(['gfortran', '-fdefault-real-8', '-O3', _run_file_])
    if _run_file_ == None:
        _run_file_ = glob.glob('*.f95')
        _run_file_ = _run_file_[0]
    # fortran_file = glob.glob('*.f95')
    p = call(['gfortran', '-fdefault-real-8', '-O3', _run_file_])
    p = call(['./a.out'])

    # remove executable and module files
    p = call(['rm', 'a.out'])
    for fl in glob.glob('*.mod'):
        os.remove(fl)


cwd = os.getcwd()
print(cwd)

CORE_MODULES_PATH = '/Users/brady22/Documents/Aestec/Models/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods'
# run_fortran()

fortran_template = 'Diffusion_1Specie_DNAD_template.f95'
fortran_specific_file = 'Diffusion_1Specie_DNAD_test.f95'

# Write the correct CORE_MODULES_PATH in the file - this is in lines that use "include"
sed(pattern='__CORE__MODULES_PATH__', 
    replace=CORE_MODULES_PATH, 
    source=fortran_template, dest=fortran_specific_file)

# run the fortran code
run_fortran(fortran_specific_file)
# sort the modules in the fortran file
# p = call(['SortFortranModules.py', fortran_template])
# _sorted_fortran_template = 'ReSort' + fortran_template
# run_fortran(_sorted_fortran_template)        # run fortran with values
