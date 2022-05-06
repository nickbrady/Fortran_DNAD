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
import time

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
    p = call(['gfortran', '-fdefault-real-8', '-O3', _run_file_])
    p = call(['./a.out'])

    # remove executable and module files
    p = call(['rm', 'a.out'])
    for fl in glob.glob('*.mod'):
        os.remove(fl)


cwd = os.getcwd()

fortran_template = 'Electrode_HalfCell.f95'
fortran_values = 'Electrode_HalfCell_values.f95'

# sort the modules in the fortran file
#p = call(['SortFortranModules.py', fortran_template])
#_sorted_fortran_template = 'ReSort' + fortran_template

c_rates = [1, 2, 3, 4, 5, 10, 20]#[10, 5, 4, 3, 2, 1]
for c_rate in c_rates:

    _source_tran = fortran_template #_sorted_fortran_template#'{}/{}'.format(cwd, _sorted_fortran_template)
    _temp = 'fortran_temp.f95'
    _dest_fold = cwd + '/{}'.format(c_rate)

    if not(os.path.exists(_dest_fold)):   # if folder does not exist, create it
        os.makedirs(_dest_fold)

    print(_dest_fold)

    sed('__C_RATE__', str(c_rate), _source_tran, _temp)
    shutil.move( _temp, '{}/{}'.format(_dest_fold, fortran_values) )

    os.chdir(_dest_fold)               # migrate to directory
    t0 = time.time()
    run_fortran(fortran_values)        # run fortran with values
    print(time.time() - t0, ' seconds')

    # remove temporary files
    os.chdir(cwd)
    for fl in glob.glob('*temp*'):
    	os.remove(fl)
