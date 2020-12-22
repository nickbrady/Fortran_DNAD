import os
import pandas as pd

def Fortran_Module_DataFrame(filename):
	'''
	Input  - fortran file
	Output - DataFrame
				* ModuleName
				* StartLine
				* EndLine
				* Len_UseList 	- number of modules in the use statements
				* Use_List 		- list of modules in the use statements
		- DataFrame is returned sorted by number of modules it "uses"
	'''
	# use rstrip to remove '\n' from end of each line
	# use .lower() because fortran is not case-sensitive
	with open(filename) as f:
	    Lines = [line.rstrip().lower() for line in f]

	module_df = pd.DataFrame(columns = ['ModuleName', 'StartLine', 'EndLine', 'Len_UseList', 'Use_List'])

	module_name = ''
	mod_num = 0
	in_module = False
	for line_num, line_str in enumerate(Lines):
		comment_position = line_str.find("!")		# -1 means ! not found
		# comment_position = min(line_str.find("c"), comment_position)
		if comment_position == 0: # ignore comment lines
			continue

		# consider linestring only before comment
		if comment_position != -1:
			line_str = line_str[:comment_position]

		# module needs to appear on its own (not "end module" or "module procedure")
		if "module " in line_str:
			if "end module " not in line_str and "module procedure " not in line_str:
				# get module name, start line, and initiate use_list
				module_name = line_str[line_str.index("module ") + len("module "):]
				module_start = line_num
				use_list = []
				in_module = True

		if "use " in line_str and in_module: 									# use statement
			# module names contained between "use" and ",only" separated by ","
			# remove whitespace to simplify pattern recognition
			line_str = line_str[line_str.index("use ") + len("use "):]
			line_str = line_str.replace(' ', '')			# remove whitespace
			only_position = line_str.find(",only")
			if only_position != -1:
				line_str = line_str[:only_position]

			use_mods = line_str.split(',')
			for mod in use_mods:
				if mod not in use_list:	# only add modules not already listed
					use_list.append(mod)

		# find the end of the module and add information to dataframe
		if module_name != '':
			if "end module " + module_name in line_str: 			# end of module
				module_end = line_num
				module_df.loc[mod_num] = [module_name, module_start, module_end, len(use_list), use_list]
				mod_num += 1
				in_module = False

	module_df = module_df.sort_values('Len_UseList')
	return module_df.reset_index(drop=True)



def sortFortranModules(ModuleDataFrame):
	'''
	Input  - ModuleDataFrame - sorted by the number of used modules
	Output - ModuleDataFrame sort in an order of the modules that will allow compilation in gfortran
	'''
	# sort the modules so that gfortran can compile
	# dependent modules need to appear after parent modules
	num_parent_modules = len(ModuleDataFrame.loc[ModuleDataFrame['Len_UseList'] == 0])
	parent_mods = list(ModuleDataFrame.loc[ModuleDataFrame['Len_UseList'] == 0]['ModuleName'].values)

	for i in range(num_parent_modules, len(ModuleDataFrame)):
		for j in range(i, len(ModuleDataFrame)):
			# if all the use modules are in parent_mods then the module can be compiled - move the module to row#: 1+num_parent_modules
			if set(ModuleDataFrame.loc[j]['Use_List']).issubset(set(parent_mods)):
				if j > num_parent_modules:
					# exchange rows
					temp = ModuleDataFrame.loc[num_parent_modules].copy()
					ModuleDataFrame.loc[num_parent_modules] = ModuleDataFrame.loc[j]
					ModuleDataFrame.loc[j] = temp
				parent_mods.append(ModuleDataFrame.loc[num_parent_modules]['ModuleName'])
				num_parent_modules += 1

				break

	return ModuleDataFrame





def ReWriteFortran(FortranFile, OutPutFile=None):
	'''
	Input - FortranFile(.f95)
				* file with modules not in order necessary for compilation in gfortran
	Output - SortedFortranFile(.f95)
				* file with modules in order necessary for compilation in gfortran

		FortranFile > Fortran_Module_DataFrame(FortranFile) > ModuleDataFrame
		sortFortranModules(ModuleDataFrame) > sorted_ModuleDataFrame
		WriteSortedFortranFile(sorted_ModuleDataFrame) > SortedFortranFile
	'''

	module_df = Fortran_Module_DataFrame(FortranFile)
	module_df = sortFortranModules(module_df)

	file1 = open(FortranFile, 'r')
	Lines = file1.readlines()

	# use module_df to re-order lines in fortran file
	linesToWrite = []
	for index, row in module_df.iterrows():
		linesToWrite = linesToWrite + list(range(row['StartLine'], row['EndLine'] + 1))

	if OutPutFile == None:
		OutPutFile = 'ReSort' + FortranFile
	with open(OutPutFile, 'w') as SortedFortranFile:
	    # write the modules first
		for idx in linesToWrite:
			SortedFortranFile.write(Lines[idx])

			if idx in module_df['EndLine'].values:
				SortedFortranFile.write('!*'*40 + '\n')
				SortedFortranFile.write('\n')

		# write the rest of the file
		for i in range(len(Lines)):
			if i in linesToWrite:
				continue
			SortedFortranFile.write(Lines[i])


# In[2]:
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/')

module_df = Fortran_Module_DataFrame("DNAD_BAND.f95")

module_df


# In[3]:
# sort the modules so that gfortran can compile
# dependent modules need to appear after parent modules
num_parent_modules = len(module_df.loc[module_df['Len_UseList'] == 0])
parent_mods = list(module_df.loc[module_df['Len_UseList'] == 0]['ModuleName'].values)

for i in range(num_parent_modules, len(module_df)):
	for j in range(i, len(module_df)):
		# if all the use modules are in parent_mods then the module can be compiled - move the module under the
		if set(module_df.loc[j]['Use_List']).issubset(set(parent_mods)):
			print(module_df.loc[j]['ModuleName'])
			if j != i:
				# exchange rows
				temp = module_df.loc[i].copy()
				module_df.loc[i] = module_df.loc[j]
				module_df.loc[j] = temp
			num_parent_modules += 1
			parent_mods.append(module_df.loc[j]['ModuleName'])

			break


module_df

# In[4]:

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/MakeFiles/Test/')

module_df = Fortran_Module_DataFrame("TEST.f95")
module_df

# In[5]:
module_df = Fortran_Module_DataFrame("TEST.f95")



module_df = sortFortranModules(module_df)
module_df

# In[10]:

filename = 'TEST.f95'
file1 = open(filename, 'r')
Lines = file1.readlines()

# use module_df to re-order lines in fortran file
linesToWrite = []
for index, row in module_df.iterrows():
	linesToWrite = linesToWrite + list(range(row['StartLine'], row['EndLine'] + 1))


with open('ReSort' + filename, 'w') as ReSort_file:
    # write the modules first
	for idx in linesToWrite:
		ReSort_file.write(Lines[idx])

		if idx in module_df['EndLine'].values:
			ReSort_file.write('!*'*40 + '\n')
			ReSort_file.write('\n')

	# write the rest of the file
	for i in range(len(Lines)):
		if i in linesToWrite:
			continue
		ReSort_file.write(Lines[i])

# In[11]:


os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ReSortExamples/')

ReWriteFortran('TEST_unSorted.f95')

# In[13]:

name = 'DNAD_BAND_unSorted.f95'

ReWriteFortran(name)
