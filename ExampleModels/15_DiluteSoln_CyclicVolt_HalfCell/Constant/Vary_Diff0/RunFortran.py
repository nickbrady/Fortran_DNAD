#
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

def run_fortran(_run_file_):
    p = call(['gfortran', '-fdefault-real-8', _run_file_])
    p = call(['./a.out'])

    # remove executable and module files
    p = call(['rm', 'a.out'])
    for fl in glob.glob('*.mod'):
        os.remove(fl)


cwd = os.getcwd()

fortran_template = 'DiluteSoln_CyclicVolt_HalfCell.f95'
_temp_constant = 'fortran_temp_constant_diff.f95'
_temp_diff = 'fortran_temp_diff_value.f95'

dir_diff = ''

diffusion_coef = 'Constant' #'Constant',
sed('__diff_coef__', diffusion_coef, fortran_template, _temp_constant)

diff_value = ['1e-5']#, '2e-5', '1e-6']

for d0d0 in diff_value:
    os.chdir(cwd)

    dir_diff = '{}/'.format(d0d0)

    if not(os.path.exists(dir_diff)):   # if folder does not exist, create it
        os.makedirs(dir_diff)

    sed('__d0d0__', d0d0, _temp_constant, _temp_diff)
    shutil.move( _temp_diff, dir_diff + fortran_template )
    os.chdir(dir_diff)               # migrate to directory

    # sort the modules in the fortran file
    p = call(['SortFortranModules.py', fortran_template])
    _sorted_fortran_template = 'ReSort' + fortran_template

    run_fortran(_sorted_fortran_template)        # run fortran with values

os.chdir(cwd)
for fl in glob.glob('*temp*'):
    os.remove(fl)