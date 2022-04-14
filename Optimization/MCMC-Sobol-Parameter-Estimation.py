'''
Modified April 12, 2022
Nick Brady

Christian Gould
May 7th, 2018

Python 3

Description:
	This program will estimate model parameters using Markov chain Monte Carlo method and sobol sampling.

The user must specifiy:

	Fortran_Code       		 = Name of Fortran file (.f95). Should be a folder labeled as "Model".
	Fortran_Para       		 = Name of the parameters (from your Fortran Code) you wish to estimate
	Bounds             		 = Bounds corresponding to the Fortran Parameters
	State 					 = The states you wish to test "D" "R" "C"
	Fortran_Output_File      = Your Fortran code should write the model data to a text file.
	                     	   This is the name of the text file your Fortan code writes. Columns
	                      	   headers should be the first row, followed by data.

	Output_x_column     	 = Column name header for your x-values
	Output_y_column     	 = Column name header for your y-values
	Output_state_column 	 = Column name header for your state markers. Should be in format of "R", "C", or "D"
	                     	   for Rest, Charge, or Discharge


	Exp_Files                = A list of Experimental csv files ['filename1.csv', 'filename2.csv', ...]. Should
	                           be in a folder labeled "Experiments". Column headers should be the first row
	                           followed by data.
	Exp_x_Column             = Column name header for your x-values
	Exp_y_Column             = Column name header for your y-values
	Exp_time_Column          = Column corresponding to elapsed time in hours.
	Exp_state_column         = Column name header for your state markers. Should be in format of "R", "C", or "D"
	Exp_mAhg_column          = Column which contains values with units mAh/g

	Sobol_Points             = Number of sample points you wish to evaluate over the parameter space
	Plotted_Points           = Number of points to be compared between experiment and model. Points will
							   be linearly interpolated
'''
######################################
############# USER INPUTS ############
######################################

Fortran_Code             = 'Vanadium_Just_550_Rxn_x_4.f95'
Fortran_Para             = ['dfdf','exex','nmax']
Bounds                   = [ [.01,4],[1,10],[1.2,2] ]
State                    = ['D']
Fortran_Output_File      = 'Time_Voltage.txt'

Output_x_column          = 'Time'
Output_y_column          = 'Voltage'
Output_state_column      = 'State'

Exp_Files                = ['550C_C_10.csv', '550C_C_05.csv']
Exp_y_column             = 'Volts'
Exp_x_column             = 'time h'
Exp_time_column          = 'time h'
Exp_state_column         = 'State'
Exp_mAhg_column          = 'mAh/g'

Sobol_Points             = 10
Plotted_Points           = 100

######################################
#### IMPORT NECESSARY PACKAGES #######
######################################

import numpy as np                          # Useful for numerical calculations on an array
import os                                   # Used for changing directories and open files with folders
import pandas as pd                         # Used to read csv files
import matplotlib.pyplot as plt             # Used for making plots
from sobol_seq import i4_sobol_generate     # Sobol Package used to generate Sobol Sequence
import math
from matplotlib.ticker import FormatStrFormatter
import sys
from subprocess import call                 # Allows Python to call shell script
import shutil
import random
from operator import itemgetter   			# Allows us to find groups of consecutive numbers
from itertools import groupby
import scipy.stats as stats
from scipy.optimize import minimize
from scipy.special import erf     			# error function (erf)
scipy.interpolate import interp1d

######################################
########## HELPER FUNCTIONS ##########
######################################

plt.rcParams["font.family"] = "Times New Roman"

def Cspec_Distime(File, Time_Key, mAhg_Key, State_Key):

	# Make sure all experiment files exist
	if os.path.exists(File):

		# Use Pandas to read experiment file (Column headers should be the first row)
		ExpData = pd.read_csv(File, skiprows=0)

		# Make sure all of the keys input by the user are actuallly column headers
		for key in [Time_Key, mAhg_Key, State_Key]:
			try:
				ExpData[key]
			except KeyError:
				print('** Error While Importing Experimental Data **')
				print('"{}" is not a column header in "{}"' .format(key, Exp_Filename))
				sys.exit()

		# State Column
		state_column = ExpData[State_Key]

		# Data Filtered by 'D' State
		D  = ExpData.loc[state_column == 'D']

		# mAh/g Column for Discharge
		mAhg = pd.Series.tolist(D[mAhg_Key])
		mAhg = list(map(float, mAhg))

		# Time [=] hr Column for Discharge
		time = pd.Series.tolist(D[Time_Key])
		time = list(map(float, time))

		dis_time_seconds = time[-1] * 3600
		c_specific = mAhg[-1] / time[-1]

		return(-1*c_specific, dis_time_seconds)

	else:
		print('** Error While Importing Experimental Data **')
		print('"{}" is not in {}' .format(Exp_Filename, Exp_Folder))
		sys.exit()


def Import_Experiment_Data(File, X_Key, Y_Key, State_Key, State):

	# Make sure all experiment files exist
	if os.path.exists(File):

		# Use Pandas to read experiment file (Column headers should be the first row)
		ExpData = pd.read_csv(File, skiprows=0)

		# Make sure all of the keys input by the user are actuallly column headers
		for key in [X_Key, Y_Key, State_Key]:
			try:
				ExpData[key]
			except KeyError:
				print('** Error While Importing Experimental Data **')
				print('"{}" is not a column header in "{}"' .format(key, File))
				sys.exit()

		# State Column
		state_column = ExpData[State_Key]

		# This is the data that is selected by specified state variables
		Data = ExpData[ExpData[State_Key].isin(State)]
		Data = Data.dropna(subset = [X_Key, Y_Key])

		# Convert column data to floats
		x_Data = pd.Series.tolist(Data[X_Key])
		x_Data = list(map(float, x_Data))
		y_Data = pd.Series.tolist(Data[Y_Key])
		y_Data = list(map(float, y_Data))
		return(x_Data, y_Data)

	else:
		print('** Error While Importing Experimental Data **')
		print('{} does not exist' .format(File))
		sys.exit()


def Run_Model(Input_Parameters, Model_Folder, Shell_Script):

	# PASS PARAMETERS TO SHELL SCRIPT
	cmd = ['bash', Shell_Script]
	for param in Input_Parameters:
		cmd.append(str(param))

	# RUN MODEL USING SHELL SCRIPT
	os.chdir(Model_Folder)
	p = call(cmd)


def Import_Model_Data(File, keys, State):

	# Organize Keys
	x_key     = keys[0]
	y_key     = keys[1]
	state_key = keys[2]

	# Make sure all experiment files exist
	if os.path.exists(File):

		with open(File) as temp_file:
			headers = temp_file.readlines()[0]
			headers = headers.split()

		# Use Pandas to read experiment file (Column headers should be the first row)
		SimData = pd.read_csv(File, delim_whitespace=True, names=headers, skiprows=4)

		# Make sure all of the keys input by the user are actuallly column headers
		for key in keys:
			try:
				SimData[key]
			except KeyError:
				print('** Error While Importing Model Data **')
				print('"{}" is not a column header in "{}"' .format(key, Sim_Filename))
				sys.exit()

		# This is the data that is selected by specified state variables
		Data = SimData[SimData[state_key].isin(State)]
		Data = Data.dropna(subset = [y_key, x_key])

		# Convert column data to floats
		x_Data = pd.Series.tolist(Data[x_key])
		x_Data = list(map(float, x_Data))
		y_Data = pd.Series.tolist(Data[y_key])
		y_Data = list(map(float, y_Data))

		return(x_Data, y_Data)

	else:
		print('** Error While Importing Model Data **')
		print('Model Output File: {} was not found' .format(Sim_Filename) )
		sys.exit()






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

	f_interp_sim = interp1d(SimX, SimY, kind='linear')
	SimY_matched = f_interp_sim(ExpX)

	return (ExpX, SimY_matched)


# def Linear_Interpolation(x1, y1, x2, y2, x):
# 	m = (y2-y1)/(x2-x1)
# 	return( m*(x-x1)+y1 )

# def Reduce_Data(Data1, Data2):
#
# 	SimX = Data1
# 	ExpX = Data2[0]
# 	ExpY = Data2[1]
# 	NewX = []
# 	NewY = []
#
# 	# Data must be linearly interpolated in order to have
# 	# the same x-points as xx in SimX
# 	for xx in SimX:
#
# 		# If the x-point already exists in the data we are
# 		# manipulating, then no need for linear interpolation
# 		if xx in ExpX:
# 			x = xx
# 			y = ExpY[ExpX.index(xx)]
# 			NewX.append(x)
# 			NewY.append(y)
#
# 		else:
# 			# Find the closest x-point in the data we are manipulating
# 			diff = np.array(ExpX) - xx
#
# 			try:
# 				# Find closest point to the right
# 				m = min(i for i in diff if i > 0)
# 				index_right = np.where(diff==m)[0][0]
#
# 				# Find closest point to the left
# 				m = max(i for i in diff if i < 0)
# 				index_left = np.where(diff==m)[0][0]
#
# 				# Linear interpolate the y value, then store it
# 				x = xx
# 				y = Linear_Interpolation(ExpX[index_left], ExpY[index_left], ExpX[index_right], ExpY[index_right], xx)
# 				NewX.append(x)
# 				NewY.append(y)
#
# 			# This is the exception, where the point you are trying to linear interpolate
# 			# is at the edge, and is not between two points. At this point, it will linear
# 			# interpolate by using the two closest points from the same side (right or left)
# 			except ValueError:
#
# 				diff = np.abs(diff)
# 				all_index = diff.argsort()[:]
# 				index_left = all_index[0]
# 				index_right = all_index[1]
# 				x = xx
#
# 				if ExpX[index_left] == ExpX[index_right]:
# 					index_right = all_index[2]
#
# 				y = Linear_Interpolation(ExpX[index_left], ExpY[index_left], ExpX[index_right], ExpY[index_right], xx)
# 				NewX.append(x)
# 				NewY.append(y)
# 				continue
#
# 	return(NewX, NewY)

def Least_Squares(expy, simy):
	'''Returns least squares for two data sets assuming that
	they have the same number of points and are indexed at
	the same points. You must use Reduce_Data function before
	calculating Least Squares.'''

	expy = np.array(expy)
	simy = np.array(simy)

	Diff_of_Squares = np.abs(expy - simy)**2
	LeastSquares    = np.sum(Diff_of_Squares)
	return(LeastSquares)

def CheckPaths(paths):
	'''Checks that given list of paths exists, if not it
	raises an error.'''

	for path in paths:
		if os.path.exists(path):
			pass
		else:
			os.mkdir(path)

	if os.path.exists(Fortran):
		pass
	else:
		print('** Fortran File Missing **')
		print('{} could not be found' .format(Fortran))
		sys.exit()

def MakeShellScript(Fortran_Para, Model_Folder, Fortran_Code, filename):

    fortran_name, fortran_ext = os.path.splitext(Fortran_Code)

    # Ignore this for now
    extra = ''

    with open(os.path.join(Model_Folder, filename), 'w') as shell_script:

        shell_script.write(extra + 'rm *temp*\n')
        shell_script.write(extra + 'mkdir OutputData\n')

        for i, para in enumerate(Fortran_Para):

            add = extra + 'sed s/{}/${}/g'.format(para, i+1)
            shell_script.write('{:30s}' .format(add))

            if i == 0:
                add = Fortran_Code
                shell_script.write('{:50s}' .format(add))

            else:
                add = 'fortran_temp_{}.f95' .format(chr(97+i-1))
                shell_script.write('{:50s}' .format(add))

            add = '>'
            shell_script.write('{:3s}' .format(add))

            if i == len(Fortran_Para)-1:
                shell_script.write('OutputData/{}\n' .format(fortran_name+'_WITH_VALUES'+fortran_ext))

            else:
                shell_script.write('fortran_temp_{}.f95\n' .format(chr(97+i)))

        shell_script.write('\n')
        shell_script.write('cd OutputData/\n')
        shell_script.write('gfortran {}\n' .format(fortran_name+'_WITH_VALUES'+fortran_ext) )
        shell_script.write('./a.out\n')
        shell_script.write(extra + 'rm a.out\n')
        shell_script.write(extra + 'rm *.mod\n')
        shell_script.write('cd ../\n')
        shell_script.write('sleep 0.1\n')
        shell_script.write(extra + 'rm *temp*\n')

        return(filename)

def sobol(bounds, num_points):
	'''This function returns a sobol sequence, given bounds for each
	parameter and a total number of points. Bounds for each parameter
	must be entered in a list as follows:
	bounds = [ [lwr_bnd1, upr_bnd1], [lwr_bnd2, upr_bnd2], ... ]'''

	number_of_variables = len(bounds)

	# Sobol Sequence is pseudo random from 0-1, therefore,
	# we must scale the Sobol Sequence to the bounds specified
	Starting_Sobol = i4_sobol_generate(number_of_variables, num_points)
	Modified_Sobol = Starting_Sobol

	for i in range(len(bounds)):
		add = bounds[i][0]                  # Amount to shift Sobol point in respective dimension
		diff = bounds[i][1] - bounds[i][0]  # Length to scale Sobol point in respective dimension
		Modified_Sobol[:,i] = add + Modified_Sobol[:,i] * diff

	return(Modified_Sobol)

def MCMC(mcmc_State, Exp_Files, Fortran_Para, Bounds, index):

	LVO_Sobol = pd.read_csv(os.path.join(Results_Folder, 'results.csv'))
	LVO_Sobol.head()

	Exp_Files_2 = list(Exp_Files)

	experiment = []
	for name in Exp_Files:

		header = name + ' ({})' .format(mcmc_State[0])
		c = LVO_Sobol[header]

		for i in range( len(mcmc_State) ):
			if i != 0:
				header = name + ' ({})' .format(mcmc_State[i])
				c = c + LVO_Sobol[header]
		experiment.append(c)

	total = experiment[0]
	for i in range( len(experiment) ):
		if i != 0:
			total = total + experiment[i]
	experiment.append(total)

	Exp_Files_2.append('Total')

	n_col = len(Fortran_Para)
	m_row = len(experiment)

	Histograms_Sobol = plt.figure(1, figsize=(8*n_col,5*m_row), dpi=300)

	plot_count = 1

	return_data = []

	for ex, exp_name in zip(experiment, Exp_Files_2):

		ex = ex.iloc[:index]
		ex = ex.dropna()

		## MONTE CARLO ALGORITHM
		sum_square_data = ex * 100

		len_data = len(sum_square_data)

		ran_prev = -1

		accepted_indices = []

		bootstrap_multiplier = 100

		# randomly sort the sobol sampling data
		my_randoms = random.sample(range(len_data), len_data)

		for n in range(0,len_data*bootstrap_multiplier):

			ran = random.choice(list(ex.index.values))

			if ran_prev == -1:
				accepted_indices.append(ran)
				ran_prev = ran

			else:
				# alpha is the ratio of the R_squared: if alpha > 1 then the current R_squared is better than the previous
		#             alpha = sim['R_squared'][ran_prev]/sim['R_squared'][ran]
				chi      = sum_square_data[ran]
				chi_prev = sum_square_data[ran_prev]
				alpha = np.exp( -chi/(2.) + chi_prev/(2.) )

				u = random.uniform(0,1)

				# acceptance criterion
				if alpha >= u:

					# add accepted indices to
					accepted_indices.append(ran)
					# only if the acceptance criterion is met do we update ran_prev
					ran_prev = ran

		#### PLOTTING ####
		print_data = []

		for pr, bnds in zip(Fortran_Para, Bounds):

			ax = Histograms_Sobol.add_subplot(m_row,n_col,plot_count)

			param = LVO_Sobol[pr]

			x = param[accepted_indices].values

			if len(x) == 0:
				print('ERROR')
				print('Experiment: {}, Parameter: {}' .format(exp_name, pr))
				print('With the current state values selected: {}' .format(State))
				print('Not enough of the sobol points were successfully run.')
				print('Please refer to the csv results file in the Results folder.')
				print()
				pass

			else:
				weights = np.ones_like(x)/float(len(x))

				ax.hist(x, normed=True, bins=50, weights=weights, color = 'black')

				mu, std = stats.norm.fit(x)

				# Printout of:
				# Experiment Name, Parameter Name, Mean, Standard Deviation
				print('Exp: {}, Par: {}, Mean: {:.3f}, STD: {:.3f}' .format(exp_name, pr, mu, std))

				print_data.append([mu, std])

				xmin, xmax = ax.get_xlim()
				x = np.linspace(xmin, xmax, 1000)
				p = stats.norm.pdf(x, mu, std)
				ax.plot(x, p, 'r', linewidth=2)

				ax.set_xlim(bnds)
				ax.set_xlabel(pr)
				ax.set_ylabel(exp_name)

				title = "Fit results: $\mu$ = %.3f,  $\sigma$ = %.3f" % (mu, std)
				ax.set_title(title)

				plot_count += 1

		return_data.append([index, print_data])

	state_str = ''
	for ss in mcmc_State:
		state_str = state_str + ss

	os.chdir(Results_Folder)
	plt.savefig('MCMC Plots - '+state_str)

	return(return_data)


######################################
##### ORGANIZE FILE/FOLDER PATHS #####
######################################
cwd                 = os.getcwd()
Model_Folder        = os.path.join(cwd, 'Model')
Output_Folder       = os.path.join(Model_Folder, 'OutputData')
Fortran_Output      = os.path.join(Output_Folder, Fortran_Output_File)
Fortran             = os.path.join(Model_Folder, Fortran_Code)
Exp_Folder          = os.path.join(cwd, 'Experiments')
Results_Folder      = os.path.join(cwd, 'Results')
Plots_Folder        = os.path.join(Results_Folder, 'Plots')

# Make sure Model and Experiment Folder Exist
CheckPaths([Model_Folder, Exp_Folder, Results_Folder, Plots_Folder])


######################################
####### CHECK PARAMETERS/BOUNDS ######
######################################

# Check bounds type
if len(Bounds) == len(Fortran_Para):
	pass
else:
	print('** Input Error **')
	print('Length of Bounds({}) and Parameter Names({}) Do Not Match' .format(len(Bounds), len(Fortran_Para)))
	print('Bounds:          ', Bounds)
	print('Parameter Names: ', Fortran_Para)
	sys.exit()


######################################
####### MAKE DICT OF EXP FILES #######
######################################

Experiments = {}
Experiments_Folder = os.path.join(cwd, 'Experiments')

for file in Exp_Files:
	file_path = os.path.join(Experiments_Folder, file)
	for s in State:
		Dict_Key = file+' ({})'.format(s)
		Experiments[Dict_Key] = Import_Experiment_Data(file_path, Exp_x_column, Exp_y_column, Exp_state_column, [s])


######################################
######## SET UP RESULTS FILE #########
######################################

# Generate Sobol Sequence in order to sample the area
Complete_List_of_Parameters = sobol(Bounds, Sobol_Points)

# Read the filenames in the Files Folder.
# Use these names as the index for the results dataframe.
index_list = []
for i in range(Sobol_Points):
	i += 1
	index_list.append(i)
results = pd.DataFrame({'Sobol Point':index_list})
results = results.set_index('Sobol Point')

# Use the experiment filenames and states selected to
# create appropriate columns
for p in Fortran_Para:
	results[p] = np.nan

# Use the experiment filenames and states selected to
# create appropriate columns for results csv file
for exp in Exp_Files:
	for s in State:
		results[exp+' ({})'.format(s)] = np.nan
		results[exp+' ({}) # Points'.format(s)] = np.nan

# Fill in each sobol point in the results csv file
j = 0
for pt in Complete_List_of_Parameters:
	j += 1
	for parameter, p in zip(pt, Fortran_Para):
			results[p][j] = parameter

# Headers for pandas file
final_headers = []
for p in Fortran_Para:
	final_headers.append(p)
for exp in Exp_Files:
	for s in State:
			final_headers.append(exp+' ({})'.format(s))
			final_headers.append('Points')

######################################
########## RUN MAIN PROGRAM ##########
######################################

Sim_Keys = [Output_x_column, Output_y_column, Output_state_column]
Shell_Script = MakeShellScript(Fortran_Para + ['cspec', 'distime'], Model_Folder, Fortran_Code, 'shell_script.sh')

# Sobol Index starting at 0
index = 0

# Iterate through the list of sobol parameters that have been premade
for Input_Parameters in Complete_List_of_Parameters:

	# Create index for each sobol point
	index += 1

	# Edit Sobol_Start if you wish to start at a later value instead of 0
	Sobol_Start = 0
	if index >= Sobol_Start:

		# Iterate through the experiments
		for exp_name in Exp_Files:

			# File path to experiment file
			path = os.path.join(Experiments_Folder, exp_name)

			# These values are calculated from the experiment file
			# and used as inputs to run the Fortran program
			cspec, distime = Cspec_Distime(path, Exp_time_column, Exp_mAhg_column, Exp_state_column)

			# Run the Fortran code with input parameters
			Run_Model(Input_Parameters + [cspec] + [distime], Model_Folder, Shell_Script)

			# At this point, a text file has been made which contains the data
			# from the Fortran Model that is marked by state, "D" "R" "C"
			# In this for loop you will calculate the error for each state
			for k, s in enumerate(State):

				# Import the Fortran model data
				SimX, SimY = Import_Model_Data(Fortran_Output, Sim_Keys, [s])

				# If the data has less than 5 points, ignore it because it is not useful
				if len(SimX) < 5:
					SimX = [0]
					SimY = [0]
					ExpX = [0]
					ExpY = [0]
					Reduced_SimX = [0]
					Reduced_SimY = [0]
					Reduced_ExpX = [0]
					Reduced_ExpY = [0]
					error = 'No Value'
					continue

				else:
					# Open corresponding experiment data that was organized in
					# a dictionary earlier
					Dict_Key = exp_name + ' ({})'.format(s)
					ExpX, ExpY = Experiments[Dict_Key]

					# Experiment and model data do not have the same points
					# The x-range of the experiment and model data will be
					# determined by the model data. At this point, the x-range of
					# of the model data will divided into a certain number of points
					# specified by the user.
					x_range = np.linspace(SimX[0], SimX[-1], Plotted_Points)

					# Use linear interpolation to reduce the experiment and model data
					# down to points that have the x-range values specified above
					Reduced_ExpX, Reduced_ExpY = Reduce_Data( x_range , (ExpX, ExpY) )    # Experiment data fit to x-range steps
					Reduced_SimX, Reduced_SimY = Reduce_Data( x_range , (SimX, SimY) )    # Model data fit to x-range steps

					# Calculate the least squares error value for the current state in the loop
					error = Least_Squares(Reduced_ExpY, Reduced_SimY)

				# Save Error Calculation to pandas table
				results[exp_name+' ({})'.format(s)][index]          = error
				results[exp_name+' ({}) # Points'.format(s)][index] = len(Reduced_SimX)

				os.chdir(Results_Folder)
				results.to_csv('results.csv', header=final_headers)

# Run MCMC
MCMC(State, Exp_Files, Fortran_Para, Bounds, Sobol_Points)
