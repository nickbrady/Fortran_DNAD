'''
Written by Nicholas Brady
March 16, 2021
'''

# In[0]:
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().run_line_magic('precision', '16')
import numpy as np                          # Useful for numerical calculations on an array
import os                                   # Used for changing directories and open files with folders
import pandas as pd                         # Used to read csv files
import matplotlib
import matplotlib.pyplot as plt             # Used for making plots
from matplotlib.ticker import FormatStrFormatter
import sys

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

# print(os.getcwd())
cw = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/'
os.chdir(cw)

# In[1]:
# ------------------------------------------------------------------------------
# Simple Diffusion (No Rxn) - 1 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('Diffusion')

# In[2]:
# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

for time in [0, 0.1, 1, 3, 10]:
    plt.plot(data['Position'].where(data['Time'] == time),
             data['Conc'].where(data['Time'] == time))

plt.ylim(0, 1.1)
plt.title('Concentration Profiles vs Time')
plt.ylabel('Concentration (mol/L)')
plt.xlabel('Position (μm)')






# In[3]:
# ------------------------------------------------------------------------------
# Diffusion Reaction - 1 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('Diffusion_Reaction')

# In[4]:
# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

for time in [0, 0.1, 1, 3, 10]:
    plt.plot(data['Position'].where(data['Time'] == time),
             data['Conc'].where(data['Time'] == time))

plt.ylim(0, 1.1)
plt.title('Concentration Profiles vs Time')
plt.ylabel('Concentration (mol/L)')
plt.xlabel('Position (μm)')


# In[5]:
# ------------------------------------------------------------------------------
# Diffusion Reaction - 2 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('TwoSpecie_Diff_Rxn')


# In[6]:
# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)

for time in [0, 0.1, 1, 3, 10]:
    ax[1].plot(data['Position'].where(data['Time'] == time),
               data['Conc_0'].where(data['Time'] == time))
    ax[2].plot(data['Position'].where(data['Time'] == time),
               data['Conc_x'].where(data['Time'] == time))

ax[1].set_ylim(0, 1.1)
ax[1].set_title('c0 vs Time')
ax[1].set_ylabel('Concentration (mol/L)')
ax[1].set_xlabel('Position (μm)')

ax[2].set_ylim(0, 1.1)
ax[2].set_title('c_x vs Time')
ax[2].set_ylabel('Concentration (mol/L)')
ax[2].set_xlabel('Position (μm)')

# In[5]:
# ------------------------------------------------------------------------------
# Electrode Model - Half Cell
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('Electrode_HallCell')


# In[6]:
# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
rows, columns = [2, 2]
ax2, fig2 = axes(2, rows, columns)

df_slice = data[(data['Position'] == data['Position'].max()) & (data['Voltage'] <= 3.0)]

NJ = len(np.unique(data['Position']))

ax[1].plot(df_slice['Equivalence'], df_slice['Voltage'])

ax[1].set_ylabel('Voltage (Volts)')
ax[1].set_xlabel('Equivalence $\mathregular{Li_xFe_3O_4}$')


df_slice_D = data[(data['State'] == 'D')]
equivs = np.unique(df_slice_D['Equivalence'])


for equiv in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
    min_ind = np.where(abs(equivs - equiv) == abs(equivs - equiv).min())
    equiv = equivs[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Equivalence'] == equiv)]
    ax2[1].plot(df_slice_D['Position'], df_slice_D['Soln_Conc'])
    ax2[3].plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])
    ax2[2].plot(df_slice_D['Position'], df_slice_D['Voltage'])
    ax2[4].plot(df_slice_D['Position'], df_slice_D['Solution_Pot'])


ax2[1].set_title('Solution Concentration')
ax2[1].set_xlabel('Position (μm)')
ax2[1].set_ylabel('Concentration (mol/L)')

ax2[3].set_title('Solid-State Concentration')
ax2[3].set_xlabel('Position (μm)')
ax2[3].set_ylabel('Equivalence $\mathregular{Li_xFe_3O_4}$')

ax2[2].set_title('Solid-State Potential')
ax2[2].set_xlabel('Position (μm)')
ax2[2].set_ylabel('Potential (Volts)')

ax2[4].set_title('Solution Potential')
ax2[4].set_xlabel('Position (μm)')
ax2[4].set_ylabel('Potential (Volts)')

fig2.tight_layout()
