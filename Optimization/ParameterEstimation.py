get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('precision', '16')
import numpy as np                          # Useful for numerical calculations on an array
import os                                   # Used for changing directories and open files with folders
import pandas as pd                         # Used to read csv files
import matplotlib
import matplotlib.pyplot as plt             # Used for making plots
from matplotlib.ticker import FormatStrFormatter
import sys
from subprocess import call                 # Allows Python to call shell script
import shutil  								# allows python to execute terminal commands (such as rm mv cp)
import glob 								# module used to find files in a directory (used to remove
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy import interpolate
from sobol_seq import i4_sobol_generate     # Sobol Package used to generate Sobol Sequence
import re

import tempfile as tp

def plot_parameters():
    matplotlib.rcParams['lines.linewidth'] = 3    # change the default line thickness for plots to 4
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 15}            # change the default font size and type for axes labels and ticks
    plt.rc('font', **font)

    csfont = {'fontname':'Serif'}
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'

matplotlib.rc('xtick', top=True, bottom=True, direction='in')
matplotlib.rc('ytick', left=True, right=True, direction='in')
plt.rc("axes.spines", top=True, right=True)
matplotlib.rc('axes', edgecolor='k')

plot_parameters()


def cell_width(Percentage):
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:" + str(Percentage) + "% !important; }</style>")) # set the jupyter notebook cell width to xx%

cell_width(90)


def axes(number, rows, columns):
    fig = plt.figure(number, figsize=(6*columns, 5*rows), dpi=80)
    ax = []
    ax.append([])
    for i in range(1, rows*columns+1):
        ax.append(fig.add_subplot(rows, columns, i))

    return ax, fig

def sobol(bounds, num_points):
	'''This function returns a sobol sequence, given bounds for each
	parameter and a total number of points. Bounds for each parameter
	must be entered in a list as follows:
	bounds = [ [lwr_bnd1, upr_bnd1], [lwr_bnd2, upr_bnd2], ... ]'''

	number_of_variables = len(bounds)

	# Sobol Sequence is pseudo random from 0-1, therefore,
	# we must scale the Sobol Sequence to the bounds specified
	Sobol_Points = i4_sobol_generate(number_of_variables, num_points)

	for i in range(len(bounds)):
		add = bounds[i][0]                  # Amount to shift Sobol point in respective dimension
		diff = bounds[i][1] - bounds[i][0]  # Length to scale Sobol point in respective dimension
		Sobol_Points[:,i] = add + Sobol_Points[:,i] * diff       # Modify Sobol point in respective dimension

	return(Sobol_Points)

# In[1]:
'''
	True function: 2x² + x + 3
'''

np.random.seed(1)

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

function_master = lambda x: 2 * x**2 + 1 * x + 3

x = np.linspace(0, 3, 100)
noise = np.random.normal(0, 1, len(x))

ax[1].plot(x, function_master(x), 'r-', linewidth=2)
ax[1].plot(x, function_master(x) + noise, 'ko', fillstyle='none')

# In[2]:
# use minimize to find coefficients of function_master
def error(vector_of_coefficients, xdata, ydata):

    estimated_values = np.polyval(vector_of_coefficients, xdata)

    return sum( (ydata - estimated_values)**2 )

coeffs_0 = [0, 0, 0]
res = minimize(error, coeffs_0, args = (x_data, y_true), options={'disp': True}) # method=BFGS
print(res.x)

res = minimize(error, coeffs_0, args = (x_data, y_true), method='nelder-mead',
               options={'disp': True})

# In[2]:
np.random.seed(1)

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
x_data = np.linspace(0, 3, 200)
y_true = function_master(x)
noise = np.random.normal(0, 1, len(x))
y_real = y_true + noise
ax[1].plot(x_data, y_true, 'r-', linewidth=5)
ax[1].plot(x_data, y_real, 'ko', fillstyle='none')
#
# use minimize to find coefficients of function_master
def error(vector_of_coefficients, xdata, ydata):

    estimated_values = np.polyval(vector_of_coefficients, xdata)

    return sum( (ydata - estimated_values)**2 )
#
print("minimum error:", sum(noise**2))
print("error with true coefficients:", error([2, 1, 3], x_data, y_real))

coeffs_0 = np.random.uniform(size=6)
res = minimize(error, coeffs_0, args = (x_data, y_true), options={'disp': True})
print("Even using a 6th-order polynomial, the derived coefficients are still very close to the true values:")
print(res.x)
# print(res)
print()


for p in range(1,7):
	coeffs_0 = np.random.uniform(size=p)
	res = minimize(error, coeffs_0, args = (x_data, y_real))
	print(res.message)
	print('Function Value:', res.fun, '\t Function Evals:', res.nfev)
	print(res.x)
	ax[1].plot(x, np.polyval(res.x, x), '--', linewidth=2)
	ax[2].semilogy(p, res.fun, 's')
	print()

# In[3]:
'''
	Nelder-Mead method produces very similar results to the BFGS method
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax[1].plot(x_data, y_true, 'r-', linewidth=5)
ax[1].plot(x_data, y_real, 'ko', fillstyle='none')

for p in range(1,7):
	coeffs_0 = np.random.uniform(size=p)
	res = minimize(error, coeffs_0, args = (x_data, y_real), method='nelder-mead')
	print(res.message)
	print('Function Value:', res.fun, '\t Function Evals:', res.nfev)
	print(res.x)
	ax[1].plot(x, np.polyval(res.x, x), '--', linewidth=2)
	ax[2].semilogy(p, res.fun, 's')
	print()

# In[4]:
# np.random.seed(1)
#
# rows, columns = [1, 1]
# ax, fig = axes(1, rows, columns)
# x = np.linspace(0, 3, 200)
# noise = np.random.normal(0, 1, len(x))
# ax[1].plot(x, function_master(x), 'r-', linewidth=2)
# ax[1].plot(x, function_master(x) + noise, 'ko', fillstyle='none')
#
#
#
# Bounds       = [ [0, 10], [-10, 10], [-10, 10] ]
# num_points = 1000
# Complete_List_of_Parameters = sobol(Bounds, num_points)
# error_vals = np.zeros(len(Complete_List_of_Parameters))
#
# for o, _params_ in enumerate(Complete_List_of_Parameters):
# 	error_vals[o] = error(_params_)
#
# df = pd.DataFrame(Complete_List_of_Parameters)
# df['Error'] = error_vals
#
# print(df.sort_values('Error'))
# print(df['Error'].idxmin())
# df.iloc[742].values[:3]
# res = minimize(error, df.iloc[742].values[:3], options={'disp': True})
# print(res.x)
# ax[1].plot(x, np.polyval(res.x, x), '-', color = 'orange', linewidth=2)
#
# res = minimize(error, df.iloc[742].values[:3], method='nelder-mead', options={'disp': True})
# print(res.x)
# ax[1].plot(x, np.polyval(res.x, x), '-', color = 'green', linewidth=2)
# print(error(res.x))
#
# res = minimize(error, res.x, options={'disp': True})
# print(res.x)
# ax[1].plot(x, np.polyval(res.x, x), '-', color = 'orange', linewidth=2)
#
#
# print(error([2,1,3]))
#
# res = minimize(error, [2,1,3], options={'disp': True})
# print(res.x)

# In[10]:
'''
	Estimate D_LiPF6 in EC:EMC
'''
exp_data_directory = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_LiPF6_in_EC_EMC/With_Current/'
sim_data_directory = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D/With_Current/'

exp_Data_File = 'Time_Voltage_Position.txt'
exp_data = pd.read_table(exp_data_directory + exp_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(exp_data.keys())
max_values = exp_data.max()
exp_Delta_Phi_2 = exp_data['Φ_2'][exp_data['Position'] == 0].values \
			  	  - exp_data['Φ_2'][exp_data['Position'] == max_values.Position].values
exp_time = np.unique(exp_data['Time'].values)

plt.plot(exp_time, exp_Delta_Phi_2, 'k--')


sim_Data_File = 'Time_Voltage_Position.txt'
sim_data = pd.read_table(sim_data_directory + sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(sim_data.keys())
max_values = sim_data.max()
sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
			  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
sim_time = np.unique(sim_data['Time'].values)

plt.plot(sim_time, sim_Delta_Phi_2, 'r-')

# In[11]:
def Reduce_Data(Simulated_Data, Experimental_Data):
	'''
		The points from the simulated and real data do not exactly align.
		For example, experimental data could be taken at 1, 1.5, 1.63, ... (seconds)
		and the simulated data could just be at 0, 1, 2, 3, ...
		We align our points by creating a function to linearly interpolate our
		simulated data - f_interp_sim. Then f_interp_sim(ExpX) produces the simulated
		data values at the exact x-values in the experimental data
	'''

	SimX = Simulated_Data[0]
	SimY = Simulated_Data[1]

	ExpX = Experimental_Data[0]
	ExpY = Experimental_Data[1]

	f_interp_sim = interp1d(SimX, SimY, kind='linear', fill_value="extrapolate")
	SimY_matched = f_interp_sim(ExpX)

	return (ExpX, SimY_matched)

Sim_data_Reduced = Reduce_Data((sim_time, Delta_Phi_2), (exp_time, exp_Delta_Phi_2))


plt.plot(Sim_data_Reduced[0], Sim_data_Reduced[1], 'r-')
print(sum((Sim_data_Reduced[1] - exp_Delta_Phi_2)**2))

# In[14]:
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

# In[13]:
wd = "/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D/With_Current/"
os.chdir(wd)

source = "LiPF6_in_EC_EMC.f95"
destination = "LiPF6_in_EC_EMC_val.f95"
sed('__p1__', str(2.0e3), wd+source, dest=wd+destination)

# In[14]:


p = call(['gfortran', '-fdefault-real-8', '-O3', destination])
p = call(['./a.out'])

p = call(['rm', 'a.out'])
for fl in glob.glob('*.mod'):
	os.remove(fl)

# In[15]:
sim_Data_File = 'Time_Voltage_Position.txt'
sim_data = pd.read_table(wd+sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
max_values = exp_data.max()
sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
			  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
sim_time = np.unique(sim_data['Time'].values)

Sim_data_Reduced = Reduce_Data((sim_time, sim_Delta_Phi_2), (exp_time, exp_Delta_Phi_2))
error = sum((Sim_data_Reduced[1] - exp_Delta_Phi_2)**2)
print(error)

plt.plot(sim_time, sim_Delta_Phi_2, 'r-')
# In[16]:
wd = "/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D/With_Current/"
os.chdir(wd)

source = "LiPF6_in_EC_EMC.f95"
destination = "LiPF6_in_EC_EMC_val.f95"

def parameter_estimation_error(params):
	print(params)
	_source = wd+source
	_destination = wd+destination
	sed('__p1__', str(params[0]), _source, dest=_destination)
	p = call(['gfortran', '-fdefault-real-8', '-O3', _destination])
	p = call(['./a.out'])

	p = call(['rm', 'a.out'])
	for fl in glob.glob('*.mod'):
		os.remove(fl)

	sim_Data_File = 'Time_Voltage_Position.txt'
	sim_data = pd.read_table(wd+sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
	max_values = sim_data.max()
	sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
				  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
	sim_time = np.unique(sim_data['Time'].values)

	Sim_data_Reduced = Reduce_Data((sim_time, sim_Delta_Phi_2), (exp_time, exp_Delta_Phi_2))
	error = sum((Sim_data_Reduced[1] - exp_Delta_Phi_2)**2)
	print(error)

	return error

params_0 = 1e3
res = minimize(parameter_estimation_error, params_0, method='nelder-mead', options={'disp': True, 'maxiter': 100})

# In[20]:
'''
	Add noise to the experimental data
'''

exp_Data_File = 'Time_Voltage_Position.txt'
exp_data = pd.read_table(exp_data_directory + exp_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(exp_data.keys())
max_values = exp_data.max()
exp_Delta_Phi_2 = exp_data['Φ_2'][exp_data['Position'] == 0].values \
			  	  - exp_data['Φ_2'][exp_data['Position'] == max_values.Position].values
noise = np.random.normal(0, 1e-2, len(exp_Delta_Phi_2))

exp_time = np.unique(exp_data['Time'].values)

plt.plot(exp_time, exp_Delta_Phi_2, 'k-', linewidth=2)
plt.plot(exp_time, exp_Delta_Phi_2 + noise, 'ko', fillstyle='none')


sim_Data_File = 'Time_Voltage_Position.txt'
sim_data = pd.read_table(sim_data_directory + sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(exp_data.keys())
max_values = exp_data.max()
sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
			  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
sim_time = np.unique(sim_data['Time'].values)

plt.plot(sim_time, sim_Delta_Phi_2, 'r-')

# In[21]:
'''
	Search for 2 parameters
'''
exp_Delta_Phi_2 = exp_data['Φ_2'][exp_data['Position'] == 0].values \
			  	  - exp_data['Φ_2'][exp_data['Position'] == max_values.Position].values
noise = np.random.normal(0, 1e-2, len(exp_Delta_Phi_2))
exp_Delta_Phi_2 += noise

params_0 = 1e3
res = minimize(parameter_estimation_error, params_0, method='nelder-mead', options={'disp': True, 'maxiter': 100})

# In[22]:
exp_Delta_Phi_2 = exp_data['Φ_2'][exp_data['Position'] == 0].values \
			  	  - exp_data['Φ_2'][exp_data['Position'] == max_values.Position].values
noise = np.random.normal(0, 1e-2, len(exp_Delta_Phi_2))
exp_Delta_Phi_2 += noise

def parameter_estimation_error(params):
	print(params)
	_source = wd+source
	_destination = wd+destination

	shutil.copyfile(_source, _destination)

	parameter_strings = ['__p1__', '__p2__']
	for (p_string, p_val) in zip(parameter_strings, params):
		sed(p_string, str(p_val), _destination)

	p = call(['gfortran', '-fdefault-real-8', '-O3', destination])
	p = call(['./a.out'])

	p = call(['rm', 'a.out'])
	for fl in glob.glob('*.mod'):
		os.remove(fl)

	sim_Data_File = 'Time_Voltage_Position.txt'
	sim_data = pd.read_table(wd+sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
	max_values = exp_data.max()
	sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
				  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
	sim_time = np.unique(sim_data['Time'].values)

	Sim_data_Reduced = Reduce_Data((sim_time, sim_Delta_Phi_2), (exp_time, exp_Delta_Phi_2))
	error = sum((Sim_data_Reduced[1] - exp_Delta_Phi_2)**2)
	print(error)

	return error

params_0 = [1e3, 2]
res = minimize(parameter_estimation_error, params_0, method='nelder-mead',
	options={'disp': True, 'maxiter': 100, 'atol': 1e-6})

# In[22]:
x = np.array([1, 2, 4, 5])  # sort data points by increasing x value
y = np.array([2, 1, 4, 3])
arr = np.arange(np.amin(x), np.amax(x), 0.01)
s = interpolate.CubicSpline(x, y)

fig, ax = plt.subplots(1, 1)
ax.plot(x, y, 'bo', label='Data Point')
ax.plot(arr, s(arr), 'k-', label='Cubic Spline', lw=1, zorder=100)

fortran_code = ''
for i in range(x.shape[0] - 1):
	segment_x = np.linspace(x[i], x[i + 1], 100)
	# A (4, 100) array, where the rows contain (x-x[i])**3, (x-x[i])**2 etc.
	exp_x = (segment_x - x[i])[None, :] ** np.arange(4)[::-1, None]
	# Sum over the rows of exp_x weighted by coefficients in the ith column of s.c
	segment_y = s.c[:, i].dot(exp_x)
	ax.plot(segment_x, segment_y, label='Segment {}'.format(i), ls='--', lw=4)

	eqn = ''
	for k in range(s.c.shape[0]):
		# c[k, i] is a coefficient for (x-x[i])**(3-k) on the segment between x[i] and x[i+1]
		if np.sign(s.c[k,i]) != -1.0:
			eqn += '+'
		if (3 - k) == 0:
			eqn += ' {}'.format(s.c[k,i])
		elif (3-k) == 1:
			eqn += ' {}*(x-{}) '.format(s.c[k,i], x[i])
		else:
			eqn += ' {}*(x-{})**({}) '.format(s.c[k,i], x[i], 3-k)

	if i == 0:
		fortran_code += 'if '
	else:
		fortran_code += '\nelse if '
#
	fortran_code += '(x <= {}) then \n \t val = {}'.format(x[i+1], eqn)

print(fortran_code)

ax.legend()
plt.show()

# In[23]:
c_LiPF6 = np.linspace(0,3e-3)
Temp = 298

p1 = 1.01e+03
p2 = 1.01e+00
p3 = -1.56e+03
p4 = -4.87e+02

c_hat = c_LiPF6 * 1000
T = Temp

f_D_salt = lambda c: p1*np.exp(p2*c) * np.exp(p3/T) * np.exp(p4/T*c) * 1e-6
D_salt_true =  f_D_salt(c_hat)

plt.plot(c_hat, D_salt_true, 'ko', fillstyle='none')

c = np.array((0, 1.5, 3.0))
D_salt =  p1*np.exp(p2*c) * np.exp(p3/T) * np.exp(p4/T*c) * 1e-6

plt.plot(c, D_salt, 'rs')


s = interpolate.CubicSpline(c, D_salt)

plt.plot(c_LiPF6*1e3, s(c_LiPF6*1e3), 'k-', label='Cubic Spline', lw=1, zorder=100)

# In[24]:
# have program pick value of red points to minimize error
# def spline_error(D_values):
D_values = np.array((4, 1, .1)) * 1e-6
D_spline = interpolate.CubicSpline(c, D_values)

c_hat = np.linspace(min(c), max(c))

error = sum((f_D_salt(c_hat) - D_spline(c_hat))**2)
print(error)

plt.plot(c_hat, f_D_salt(c_hat))
plt.plot(c_hat, D_spline(c_hat))

# In[25]:
def error(D_values, c):

	D_spline = interpolate.CubicSpline(c, D_values)

	c_hat = np.linspace(min(c), max(c))

	error = sum((f_D_salt(c_hat) - D_spline(c_hat))**2)

	return error


D_values = np.array((10, 2, 1)) * 1e-6
res = minimize(error, D_values, args = (c,), options={'disp': True}) # method=BFGS
print(res.x)

plt.plot(c_hat, D_salt_true, 'k--', fillstyle='none')
plt.plot(c, res.x, 'rs', markersize=7)
D_spline_final = interpolate.CubicSpline(c, res.x)
plt.plot(c_hat, D_spline_final(c_hat), 'r-', linewidth=2)


# In[26]:
def write_Fortran_Spline(spline_object):

	x 		  = spline_object.x
	coeffs    = spline_object.c
	max_power = coeffs.shape[0] - 1

	precision = 12

	fortran_code = ''
	for i in range(x.shape[0] - 1):

		eqn = ''
		for k in range(max_power+1):
			# c[k, i] is a coefficient for (x-x[i])**(3-k) on the segment between x[i] and x[i+1]
			if (np.sign(coeffs[k,i]) != -1.0) and (k != 0) :
				eqn += '+'
			if (max_power - k) == 0:
				eqn += ' {:.{prec}}'.format(coeffs[k,i], prec=precision)
			elif (max_power - k) == 1:
				eqn += ' {:.{prec}}*(c-{}) '.format(coeffs[k,i], x[i], prec=precision)
			else:
				eqn += ' {:.{prec}}*(c-{})**{} '.format(coeffs[k,i], x[i], max_power - k, prec=precision)

		if i != x.shape[0] - 2:
			if i == 0:
				fortran_code += 'if '
			else:
				fortran_code += '\nelse if '
			fortran_code += '(c <= {}) then'.format(x[i+1])
		else:
			fortran_code += '\nelse '

		fortran_code += '\n\t D_salt = {}'.format(eqn)

	fortran_code += '\nend if'

	return fortran_code

# print(write_Fortran_Spline(D_spline_final))
print(write_Fortran_Spline(D_spline_final))


wd = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D_spline/With_Current/'
_source = wd+source
_destination = wd+destination

shutil.copyfile(_source, _destination)

parameter_strings = ['__FORTRAN_SPLINE_CODE__']
# for (p_string, p_val) in zip(parameter_strings, params):
sed(parameter_strings[0], write_Fortran_Spline(D_spline_final), _destination)


plt.plot(exp_time, exp_Delta_Phi_2, 'k--')

sim_Data_File = 'Time_Voltage_Position.txt'
sim_data = pd.read_table(wd + sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(sim_data.keys())
max_values = sim_data.max()
sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
			  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
sim_time = np.unique(sim_data['Time'].values)

plt.plot(sim_time, sim_Delta_Phi_2, 'r-')

# In[26]:
wd = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D_spline/With_Current/'

def parameter_estimation_error(params):
	print(params)
	_source = wd+source
	_destination = wd+destination

	shutil.copyfile(_source, _destination)

	D_values = params
	c = np.array((0, 1.5, 3.0))
	D_spline = interpolate.CubicSpline(c, D_values)

	parameter_strings = ['__FORTRAN_SPLINE_CODE__']
	sed(parameter_strings[0], write_Fortran_Spline(D_spline), _destination)

	os.chdir(wd)
	p = call(['gfortran', '-fdefault-real-8', '-O3', _destination])
	p = call(['./a.out'])

	p = call(['rm', 'a.out'])
	for fl in glob.glob('*.mod'):
		os.remove(fl)

	sim_Data_File = 'Time_Voltage_Position.txt'
	sim_data = pd.read_table(wd+sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
	max_values = sim_data.max()
	sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
				  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
	sim_time = np.unique(sim_data['Time'].values)

	Sim_data_Reduced = Reduce_Data((sim_time, sim_Delta_Phi_2), (exp_time, exp_Delta_Phi_2))
	error = sum((Sim_data_Reduced[1] - exp_Delta_Phi_2)**2)
	print(error)

	return error

params_0 = np.array((6, 2, 1)) * 1e-6
res = minimize(parameter_estimation_error, params_0, method='nelder-mead',
	options={'disp': True, 'maxiter': 100})

# In[28]:
wd = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_LiPF6_in_EC_EMC_estimate_D_spline/With_Current/'

plt.plot(exp_time, exp_Delta_Phi_2, 'k--')

sim_Data_File = 'Time_Voltage_Position.txt'
sim_data = pd.read_table(wd + sim_Data_File, delim_whitespace=True, header=0, skiprows=[1])
print(sim_data.keys())
max_values = sim_data.max()
sim_Delta_Phi_2 = sim_data['Φ_2'][sim_data['Position'] == 0].values \
			  	  - sim_data['Φ_2'][sim_data['Position'] == max_values.Position].values
sim_time = np.unique(sim_data['Time'].values)

plt.plot(sim_time, sim_Delta_Phi_2, 'r-')

print(sum((exp_Delta_Phi_2 - sim_Delta_Phi_2)**2))
