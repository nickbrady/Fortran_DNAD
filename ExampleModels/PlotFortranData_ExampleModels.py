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

import re
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text): # allows for the more natural sorting of text (1, 2, 3, ... instead of 1, 10, 2, 3,...)
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

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

Fconst = 96485

# print(os.getcwd())
cw = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/'
os.chdir(cw)

# In[1]:
# ------------------------------------------------------------------------------
# Simple Diffusion (No Rxn) - 1 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('1_Diffusion/Rectangular/')

# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

for time in [0, 0.1, 0.2, 0.5, 1.0]:
    plt.plot(data['Position'].where(data['Time'] == time),
             data['Conc'].where(data['Time'] == time))

plt.ylim(0, 1.1)
plt.title('Concentration Profiles vs Time')
plt.ylabel('Concentration (mol/L)')
plt.xlabel('Position (μm)')






# In[2]:
# ------------------------------------------------------------------------------
# Diffusion Reaction - 1 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('2_Diffusion_Reaction')

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
# Diffusion Reaction - 2 specie
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('3_TwoSpecie_Diff_Rxn')

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

# In[4]:
# ------------------------------------------------------------------------------
# Electrode Model - Half Cell
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('4_Electrode_HalfCell')

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

# In[41]:
os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
rows, columns = [2, 2]
ax2, fig2 = axes(2, rows, columns)



df_slice = data[(data['Position'] == data['Position'].max())]

NJ = len(np.unique(data['Position']))

ax[1].plot(df_slice['Equivalence'], df_slice['Voltage'], 'black')

ax[1].set_ylabel('Voltage (Volts)')
ax[1].set_xlabel('Equivalence $\mathregular{Li_xFe_3O_4}$')


df_slice_D = data[(data['State'] == 'D')]
equivs = np.unique(df_slice_D['Equivalence'])


for equiv in [0, 0.05, 0.1, 0.15, 0.2, 0.25]:
    min_ind = np.where(abs(equivs - equiv) == abs(equivs - equiv).min())
    equiv = equivs[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Equivalence'] == equiv)]

    edge_data = df_slice_D[df_slice_D['Position'] == data['Position'].max()]
    ax[1].plot(edge_data['Equivalence'], edge_data['Voltage'], 'o', markersize = 12, markeredgecolor='k')


    ax2[1].plot(df_slice_D['Position'], df_slice_D['Soln_Conc'])
    ax2[3].plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])
    ax2[2].plot(df_slice_D['Position'], df_slice_D['Voltage'])
    ax2[4].plot(df_slice_D['Position'], df_slice_D['Solution_Pot']*1e3)


# ax2[1].set_xlabel('Position (μm)')
ax2[1].set_xticklabels([])
ax2[1].set_ylabel('Solution Concentration, $\mathregular{c_0}$ (mol/L)',
fontweight='bold')

ax2[3].set_xlabel('Position (μm)', fontweight='bold')
ax2[3].set_ylabel('Solid-State Concentration ($\mathregular{c_x / c_{x,max}}$)',
fontweight='bold')

# ax2[2].set_xlabel('Position (μm)')
ax2[2].set_xticklabels([])
ax2[2].set_ylabel('Solid-State Potential, $\mathregular{Φ_1}$ (Volts)',
fontweight='bold')

ax2[4].set_xlabel('Position (μm)', fontweight='bold')
ax2[4].set_ylabel('Solution Potential, $\mathregular{Φ_2}$ (mV)',
fontweight='bold')

fig2.tight_layout()


# In[42]:
os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/')

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

itr = iter(os.walk(os.getcwd()))
root, dirs, files = next(itr)
dirs.sort(key=natural_keys)

c_rate_colors = ['red', 'black', 'green', 'darkorange']

dirs = ['1', '2', '3', '4', '5', '10', '20']
for o, dir in enumerate(dirs):
    os.chdir(dir)
    data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
    df_edge = data[(data['Position'] == data['Position'].max())]

    ax[1].plot(df_edge['Capacity'], df_edge['Voltage'])

    os.chdir('../')

ax[1].set_ylim(ymin=2.9)

# In[42]:
def Electrode_Data_pdf(filename=None):
    f = plt.figure(figsize=(12, 5))

    ax1 = plt.subplot2grid((2, 6), (0, 0), colspan=2, rowspan=2)
    ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((2, 6), (1, 2), colspan=2)
    ax5 = plt.subplot2grid((2, 6), (1, 4), colspan=2)



    itr = iter(os.walk(os.getcwd()))
    root, dirs, files = next(itr)
    dirs.sort(key=natural_keys)

    dirs = ['10/', '5/', '3/']

    c_rate_colors = ['red', 'black', 'green']

    for o, dir in enumerate(dirs):
        os.chdir(dir)
        data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
        df_edge = data[(data['Position'] == data['Position'].max())]

        ax1.plot(df_edge['Capacity'], df_edge['Voltage'], color=c_rate_colors[o])

        os.chdir('../')

    ax1.set_ylim(ymin=2.9)

    data = pd.read_table('5/Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    ax1.set_ylabel('$\mathregular{\Phi_1 |_{x=L}}$ (V vs $\mathregular{Li/Li^+}$)', fontweight='bold')
    ax1.set_xlabel('Capacity (mAh/g)', fontweight='bold')


    df_slice_D = data[(data['State'] == 'D')]
    caps = np.unique(df_slice_D['Capacity'])


    for o, _cap_ in enumerate(range(5)):#[0, 20, 40, 60]:
        _cap_ = _cap_*30
        min_ind = np.where(abs(caps - _cap_) == abs(caps - _cap_).min())
        _cap_ = caps[min_ind][0]

        df_slice_D = data[(data['State'] == 'D') & (data['Capacity'] == _cap_)]

        edge_data = df_slice_D[df_slice_D['Position'] == data['Position'].max()]

        if o == 0:
            ax1.plot(0, 3.4, 'o', markersize = 10, markeredgecolor='k')
        else:
            ax1.plot(edge_data['Capacity'], edge_data['Voltage'], 'o', markersize = 10, markeredgecolor='k')


        ax2.plot(df_slice_D['Position'], df_slice_D['Soln_Conc'])
        ax4.plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])
        ax3.plot(df_slice_D['Position'], df_slice_D['Voltage'])
        ax5.plot(df_slice_D['Position'], df_slice_D['Solution_Pot']*1e3)

    import matplotlib.ticker as ticker
    ax1.set_ylim(ymin=3.04, ymax=3.41)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

    ax1.set_xticks(np.arange(0,201,50))
    ax1.set_xticklabels([0, '', 100, '', 200])
    ax1.text(150, 3.23, 'C/10', color='red', ha='left', va='center', fontsize=20)
    ax1.text(130, 3.07, 'C/5', color='black', ha='left', va='center', fontsize=20)
    ax1.text(67, 3.07, 'C/3', color='green', ha='right', va='center', fontsize=20)

    # electrode_scale_ticks = [0, 25, 50, 75, 100]
    # ax2[1].set_xlabel('Position (μm)')
    # ax2.set_xticks(electrode_scale_ticks)
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax2.set_xticklabels([])
    ax2.set_ylabel('$\mathregular{c_0}$ (mol/L)',
    fontweight='bold')

    ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax4.set_xlabel('Position (μm)', fontweight='bold')
    ax4.set_ylabel('$\mathregular{c_x \ / \ c_{x,max}}$',
    fontweight='bold')

    # ax2[2].set_xlabel('Position (μm)')
    ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax3.set_xticklabels([])
    ax3.set_ylabel('$\mathregular{Φ_1}$ (Volts)',
    fontweight='bold')

    ax5.yaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax5.xaxis.set_minor_locator(ticker.MultipleLocator(25))
    ax5.set_xlabel('Position (μm)', fontweight='bold')
    ax5.set_ylabel('$\mathregular{Φ_2}$ (mV)',
    fontweight='bold')

    f.text(0.025, 1, 'A)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
    f.text(0.33+0.02, 1, 'B)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
    f.text(0.66+0.02, 1, 'C)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
    f.text(0.33+0.02, 0.5+0.03, 'D)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
    f.text(0.66+0.02, 0.5+0.03, 'E)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)

    f.tight_layout()
    if filename:
        f.savefig(os.getcwd() + '/Electrode_Data.pdf', format='pdf', dpi=300, bbox_inches = "tight")


os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/')

Electrode_Data_pdf(filename = os.getcwd() + '/Electrode_Data.pdf')

# In[43]:
os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/NoDualNumbers/')

Electrode_Data_pdf()

# In[44]:
os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/')

f = plt.figure(figsize=(12, 5))

ax1 = plt.subplot2grid((2, 6), (0, 0), colspan=2, rowspan=2)
ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
ax4 = plt.subplot2grid((2, 6), (1, 2), colspan=2)
ax5 = plt.subplot2grid((2, 6), (1, 4), colspan=2)



itr = iter(os.walk(os.getcwd()))
root, dirs, files = next(itr)
dirs.sort(key=natural_keys)

dirs = ['10/', '5/', '3/']

c_rate_colors = ['red', 'black', 'green']

for o, dir in enumerate(dirs):
    os.chdir(dir)
    data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
    df_edge = data[(data['Position'] == data['Position'].max())]

    ax1.plot(df_edge['Capacity'], df_edge['Voltage'], color=c_rate_colors[o])

    os.chdir('../')

ax1.set_ylim(ymin=2.9)

data = pd.read_table('5/Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

ax1.set_ylabel('$\mathregular{\Phi_1 |_{x=L}}$ (V vs $\mathregular{Li/Li^+}$)', fontweight='bold')
ax1.set_xlabel('Capacity (mAh/g)', fontweight='bold')


df_slice_D = data[(data['State'] == 'D')]
caps = np.unique(df_slice_D['Capacity'])


for o, _cap_ in enumerate(range(5)):#[0, 20, 40, 60]:
    _cap_ = _cap_*30
    min_ind = np.where(abs(caps - _cap_) == abs(caps - _cap_).min())
    _cap_ = caps[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Capacity'] == _cap_)]

    edge_data = df_slice_D[df_slice_D['Position'] == data['Position'].max()]

    if o == 0:
        ax1.plot(0, 3.4, 'o', markersize = 10, markeredgecolor='k', zorder=100)
    else:
        ax1.plot(edge_data['Capacity'], edge_data['Voltage'], 'o', markersize = 10, markeredgecolor='k', zorder=100)


    ax2.plot(df_slice_D['Position'], df_slice_D['Soln_Conc'])
    ax4.plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])
    ax3.plot(df_slice_D['Position'], df_slice_D['Voltage'])
    ax5.plot(df_slice_D['Position'], df_slice_D['Solution_Pot']*1e3)

import matplotlib.ticker as ticker
ax1.set_ylim(ymin=3.04, ymax=3.41)
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

ax1.set_xticks(np.arange(0,201,50))
ax1.set_xticklabels([0, '', 100, '', 200])
ax1.text(150, 3.23, 'C/10', color='red', ha='left', va='center', fontsize=20)
ax1.text(130, 3.07, 'C/5', color='black', ha='left', va='center', fontsize=20)
ax1.text(67, 3.07, 'C/3', color='green', ha='right', va='center', fontsize=20)

# electrode_scale_ticks = [0, 25, 50, 75, 100]
# ax2[1].set_xlabel('Position (μm)')
# ax2.set_xticks(electrode_scale_ticks)
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.set_xticklabels([])
ax2.set_ylabel('$\mathregular{c_0}$ (mol/L)',
fontweight='bold')

ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax4.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax4.set_xlabel('Position (μm)', fontweight='bold')
ax4.set_ylabel('$\mathregular{c_x \ / \ c_{x,max}}$',
fontweight='bold')

# ax2[2].set_xlabel('Position (μm)')
ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax3.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax3.set_xticklabels([])
ax3.set_ylabel('$\mathregular{Φ_1}$ (Volts)',
fontweight='bold')

ax5.yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax5.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax5.set_xlabel('Position (μm)', fontweight='bold')
ax5.set_ylabel('$\mathregular{Φ_2}$ (mV)',
fontweight='bold')

f.text(0.025, 1, 'A)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.33+0.02, 1, 'B)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.66+0.02, 1, 'C)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.33+0.02, 0.5+0.03, 'D)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.66+0.02, 0.5+0.03, 'E)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)

f.tight_layout()

os.chdir(cw)
os.chdir('4_Electrode_HalfCell/Nernst_OCP/NoDualNumbers/')
#
# f = plt.figure(figsize=(12, 5))

# ax1 = plt.subplot2grid((2, 6), (0, 0), colspan=2, rowspan=2)
# ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
# ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
# ax4 = plt.subplot2grid((2, 6), (1, 2), colspan=2)
# ax5 = plt.subplot2grid((2, 6), (1, 4), colspan=2)



itr = iter(os.walk(os.getcwd()))
root, dirs, files = next(itr)
dirs.sort(key=natural_keys)

dirs = ['10/', '5/', '3/']

c_rate_colors = ['red', 'black', 'green']

for o, dir in enumerate(dirs):
    os.chdir(dir)
    data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
    df_edge = data[(data['Position'] == data['Position'].max())]


    c= 'w'
    ax1.plot(df_edge['Capacity'], df_edge['Voltage'], color=c, ls='--', lw=1.5)

    os.chdir('../')

ax1.set_ylim(ymin=2.9)

data = pd.read_table('5/Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

ax1.set_ylabel('$\mathregular{\Phi_1 |_{x=L}}$ (V vs $\mathregular{Li/Li^+}$)', fontweight='bold')
ax1.set_xlabel('Capacity (mAh/g)', fontweight='bold')


df_slice_D = data[(data['State'] == 'D')]
caps = np.unique(df_slice_D['Capacity'])


for o, _cap_ in enumerate(range(5)):#[0, 20, 40, 60]:
    _cap_ = _cap_*30
    min_ind = np.where(abs(caps - _cap_) == abs(caps - _cap_).min())
    _cap_ = caps[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Capacity'] == _cap_)]

    edge_data = df_slice_D[df_slice_D['Position'] == data['Position'].max()]

    # if o == 0:
    #     ax1.plot(0, 3.4, 'o', markersize = 10, markeredgecolor='k')
    # else:
    #     ax1.plot(edge_data['Capacity'], edge_data['Voltage'], 'o', markersize = 10, markeredgecolor='k')


    ax2.plot(df_slice_D['Position'], df_slice_D['Soln_Conc'], color='w', ls='--', lw=1.5)
    ax4.plot(df_slice_D['Position'], df_slice_D['Solid_Conc'], color='w', ls='--', lw=1.5)
    ax3.plot(df_slice_D['Position'], df_slice_D['Voltage'], color='w', ls='--', lw=1.5)
    ax5.plot(df_slice_D['Position'], df_slice_D['Solution_Pot']*1e3, color='w', ls='--', lw=1.5)

import matplotlib.ticker as ticker
ax1.set_ylim(ymin=3.04, ymax=3.41)
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1))

ax1.set_xticks(np.arange(0,201,50))
ax1.set_xticklabels([0, '', 100, '', 200])
ax1.text(150, 3.23, 'C/10', color='red', ha='left', va='center', fontsize=20)
ax1.text(130, 3.07, 'C/5', color='black', ha='left', va='center', fontsize=20)
ax1.text(67, 3.07, 'C/3', color='green', ha='right', va='center', fontsize=20)

# electrode_scale_ticks = [0, 25, 50, 75, 100]
# ax2[1].set_xlabel('Position (μm)')
# ax2.set_xticks(electrode_scale_ticks)
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax2.set_xticklabels([])
ax2.set_ylabel('$\mathregular{c_0}$ (mol/L)',
fontweight='bold')

ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax4.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax4.set_xlabel('Position (μm)', fontweight='bold')
ax4.set_ylabel('$\mathregular{c_x \ / \ c_{x,max}}$',
fontweight='bold')

# ax2[2].set_xlabel('Position (μm)')
ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax3.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax3.set_xticklabels([])
ax3.set_ylabel('$\mathregular{Φ_1}$ (Volts)',
fontweight='bold')

ax5.yaxis.set_minor_locator(ticker.MultipleLocator(25))
ax5.xaxis.set_minor_locator(ticker.MultipleLocator(25))
ax5.set_xlabel('Position (μm)', fontweight='bold')
ax5.set_ylabel('$\mathregular{Φ_2}$ (mV)',
fontweight='bold')

f.text(0.025, 1, 'A)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.33+0.02, 1, 'B)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.66+0.02, 1, 'C)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.33+0.02, 0.5+0.03, 'D)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)
f.text(0.66+0.02, 0.5+0.03, 'E)', va = 'top', ha = 'left', fontweight='bold', fontsize=20)

f.tight_layout()
print(os.getcwd())
f.savefig(os.getcwd() + '/Electrode_Data_Compare.pdf', format='pdf', dpi=300, bbox_inches = "tight")

# In[10]:
print((48.348 - 45.638)/45.638)
print((16.800 - 15.785)/15.785)
print((6.810 - 6.093)/6.093


# In[5]:
# ------------------------------------------------------------------------------
# Agglomerate Model
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('5_Agglomerate')

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


# In[6]:
# ------------------------------------------------------------------------------
# Crystal Model - Simple Diffusion (No concentration dependence)
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('6_Crystal/SimpleDiffusion')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)

df_slice = data[(data['Position'] == data['Position'].max()) & (data['Voltage'] <= 3.1)]

NJ = len(np.unique(data['Position']))

ax[1].plot(df_slice['Equivalence'], df_slice['Voltage'])

ax[1].set_title('Measured Potential')
ax[1].set_ylabel('Potential (Volts)')
ax[1].set_xlabel('Equivalence $\mathregular{Li_xFe_3O_4}$')


df_slice_D = data[(data['State'] == 'D')]
equivs = np.unique(df_slice_D['Equivalence'])


for equiv in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
    min_ind = np.where(abs(equivs - equiv) == abs(equivs - equiv).min())
    equiv = equivs[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Equivalence'] == equiv)]
    ax[2].plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])

ax[2].set_title('Solid-State Concentration')
ax[2].set_xlabel('Position (nm)')
ax[2].set_ylabel('Equivalence $\mathregular{Li_xFe_3O_4}$')

fig.tight_layout()

# In[7]:
# ------------------------------------------------------------------------------
# Crystal Model - Concentration Dependent Diffusion
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('6_Crystal/Conc_Depend_Diff')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)

df_slice = data[(data['Position'] == data['Position'].max()) & (data['Voltage'] <= 3.1)]

NJ = len(np.unique(data['Position']))

ax[1].plot(df_slice['Equivalence'], df_slice['Voltage'])

ax[1].set_title('Measured Potential')
ax[1].set_ylabel('Potential (Volts)')
ax[1].set_xlabel('Equivalence $\mathregular{Li_xFe_3O_4}$')


df_slice_D = data[(data['State'] == 'D')]
equivs = np.unique(df_slice_D['Equivalence'])


for equiv in np.linspace(0, 1.5, 5):
    min_ind = np.where(abs(equivs - equiv) == abs(equivs - equiv).min())
    equiv = equivs[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Equivalence'] == equiv)]
    ax[2].plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])

ax[2].set_title('Solid-State Concentration')
ax[2].set_xlabel('Position (nm)')
ax[2].set_ylabel('Equivalence $\mathregular{Li_xFe_3O_4}$')

fig.tight_layout()


# In[8]:
# ------------------------------------------------------------------------------
# Crystal Model - Phase Change
# ------------------------------------------------------------------------------
os.chdir(cw)
os.chdir('6_Crystal/PhaseChange')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

df_slice = data[(data['Position'] == data['Position'].max()) & (data['Voltage'] <= 3.1)]

NJ = len(np.unique(data['Position']))

ax[1].plot(df_slice['Equivalence'], df_slice['Voltage'])

ax[1].set_title('Measured Potential')
ax[1].set_ylabel('Potential (Volts)')
ax[1].set_xlabel('Equivalence $\mathregular{Li_xFe_3O_4}$')


df_slice_D = data[(data['State'] == 'D')]
equivs = np.unique(df_slice_D['Equivalence'])


for equiv in np.linspace(0, 2.1, 5):
    min_ind = np.where(abs(equivs - equiv) == abs(equivs - equiv).min())
    equiv = equivs[min_ind][0]

    df_slice_D = data[(data['State'] == 'D') & (data['Equivalence'] == equiv)]
    ax[2].plot(df_slice_D['Position'], df_slice_D['Solid_Conc'])
    ax[4].plot(df_slice_D['Position'], df_slice_D['Theta_beta'])

ax[2].set_title('Solid-State Concentration')
ax[2].set_xlabel('Position (nm)')
ax[2].set_ylabel('Equivalence $\mathregular{Li_xFe_3O_4}$')

ax[4].set_title('Volume Fraction β-Phase')
ax[4].set_xlabel('Position (nm)')
ax[4].set_ylabel('$\mathregular{\\theta_{\\beta}}$')

ax[3].axis('off')

fig.tight_layout()


# In[9]:
# ------------------------------------------------------------------------------
# Ionic Liquid Electrolyte
# ------------------------------------------------------------------------------


MW_LiTFSI   = 287.09
MW_EMITFSI  = 391.3
MW_Li       = 6.941
MW_TFSI     = MW_LiTFSI - MW_Li
MW_EMI      = MW_EMITFSI - MW_TFSI


os.chdir(cw)
os.chdir('7_BinarySaltILE')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

NJ = len(np.unique(data['Position']))

times_array = np.array([0.0, 0.1, 0.2, 0.5, 0.99, 1.0, 1.1, 1.2, 1.5, 1.99, 3.0, 4, 5])

for time in np.unique(data['Time']):
    if any(np.isclose(times_array - time, 0.0, atol=5e-3)):

        if time < 1.0:
            df_slice = data[(data['State'] == 'D') & (data['Time'] == time)]
            velocity = df_slice['Velocity'].values
            velocity = velocity.astype(np.float)

            ax[1].plot(df_slice['Position'], df_slice['Soln_Conc'])
            ax[3].plot(df_slice['Position'], velocity)

        if time >= 1.0:
            df_slice = data[(data['State'] == 'R') & (data['Time'] == time)]
            velocity = df_slice['Velocity'].values
            velocity = velocity.astype(np.float)

            ax[2].plot(df_slice['Position'], df_slice['Soln_Conc'])
            ax[4].plot(df_slice['Position'], velocity)


        # ax[3].plot(df_slice['Position'], df_slice['c_B'])
        # ax[4].plot(df_slice['Position'], df_slice['Density'])

        # ax[5].plot(df_slice['Position'], df_slice['Mass_Frac_A'])
        # ax[5].plot(df_slice['Position'], 1.0 - df_slice['Mass_Frac_A'])

        # ax[6].plot(df_slice['Position'], df_slice['Mol_Frac_A'])
        # ax[6].plot(df_slice['Position'], 1.0 - df_slice['Mol_Frac_A'])

ax[1].set_ylabel('$c_A$  (mol/L)')
ax[1].set_xlabel('Position (μm)')

ax[2].set_ylabel('$c_A$  (mol/L)')
ax[2].set_xlabel('Position (μm)')

ax[3].set_ylabel('$𝐯$  (nm/s)')
ax[3].set_xlabel('Position (μm)')

ax[4].set_ylabel('$𝐯$  (nm/s)')
ax[4].set_xlabel('Position (μm)')


fig.tight_layout()

# In[71]:
os.chdir(cw)
os.chdir('71_ConvectionDiffusion')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [1,1]
ax, fig = axes(1, rows, columns)

NJ = len(np.unique(data['Position']))
unique_times = np.unique(data['Time'].values)

times_array = np.array([0.0, 0.01])#, 0.1, 0.2, 0.5, 0.99, 1.0, 1.1, 1.2, 1.5, 1.99, 3.0, 4, 5])

for x in [0, 1, 3, 10, 100]:
    time = unique_times[x]
    df_slice = data[data['Time'] == time]

    ax[1].plot(df_slice['Position'], df_slice['Conc']*1e3)


ax[1].set_ylabel('$c_A$  (mol/L)')
ax[1].set_xlabel('Position (μm)')


fig.tight_layout()

# In[16]:
os.chdir(cw)
os.chdir('16_CoolingFluid_in_Pipe')

rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys(), len(np.unique(data['Time'])))

unique_times = list(np.unique(data['Time']))
plot_times = unique_times# + unique_times[20::10]
# print(plot_times)
plot_times = [plot_times[-1]]
for time in plot_times:

    df_slice = data[data['Time'] == time]

    ax[1].plot(df_slice['Position'], df_slice['Temp'], color='black')
    ax[2].plot(df_slice['Position'], df_slice['Velocity'], color='blue')
    ax2_twin = ax[2].twinx()
    ax2_twin.plot(df_slice['Position'], df_slice['Density'], 'red')
    ax[2].plot(df_slice['Position'], df_slice['Density'] * df_slice['Velocity'], 'g--', linewidth=3)


ax[1].set_xticklabels([])
ax[1].set_ylabel('Temperature (K)')
ax[2].set_ylabel('$𝐯$  (cm/s)', color='blue')
ax[2].tick_params(axis='y', colors='blue')

ax2_twin.set_ylabel('$ρ$  (g/cm³)', color='red')
ax2_twin.yaxis.label.set_color('red')
ax2_twin.spines['left'].set_color('blue')        # setting up Y-axis tick color to red
ax2_twin.spines['right'].set_color('red')        # setting up Y-axis tick color to red
ax2_twin.tick_params(axis='y', colors='red')

ax[2].text(20, 8.4, '$ρ𝐯$', ha='center', fontsize = 20, color = 'green')

ax[2].set_xlabel('Position (cm)')

fig.tight_layout()


# In[20]:
os.chdir(cw)
os.chdir('20_Simple_ConcSoln_FiniteDifference')

# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())


rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

ax[1].plot(data['Position'].where(data['Time'] == max(data['Time']) ),
         data['ω_A'].where(data['Time'] == max(data['Time'])) )

ax_t = ax[1].twinx()
ax_t.plot(data['Position'].where(data['Time'] == max(data['Time']) ),
         data['vel'].where(data['Time'] == max(data['Time'])) )

ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_A}$')
ax[1].set_xlabel('Position (cm)')

# In[21]:
os.chdir(cw)
os.chdir('21_ConcSoln_2_FiniteDifference')

# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())


rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

ax[1].plot(data['Position'].where(data['Time'] == max(data['Time']) ),
         data['ω_A'].where(data['Time'] == max(data['Time'])) )

ax_t = ax[1].twinx()
ax_t.plot(data['Position'].where(data['Time'] == max(data['Time']) ),
         data['vel'].where(data['Time'] == max(data['Time'])) )

ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_A}$')
ax[1].set_xlabel('Position (cm)')

# In[22]:
os.chdir(cw)
os.chdir('22_BinarySalt_ConcSoln')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

for time in [0, 1, 10, 30, 300]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[4].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    # if time == 0:
    #     continue
    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )



ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[2].set_title('Mass Average Velocity vs Time')
ax[2].set_ylabel('Velocity (nm/s)')
# ax[2].set_xlabel('Position (cm)')

ax[3].set_title('Flux vs Time')
ax[3].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[3].set_xlabel('Position (cm)')

ax[4].set_title('Flux vs Time')
ax[4].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[4].set_xlabel('Position (cm)')

fig.tight_layout()

# In[23]:
os.chdir(cw)
os.chdir('23_BinarySalt_ConcSoln_withCurrent')
# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

time_unique = np.unique(data['Time'])

print(data.keys())

rows, columns = [4, 3]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()


for time in [0, time_unique[1], 1, 10, 30, 300]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )

    ax[columns+2].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_+'].where(data['Time'] == time) )

    ax[columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_-'].where(data['Time'] == time) )


    ax[2*columns+1].plot(data['Position'].where(data['Time'] == time ),
       data['i_+'].where(data['Time'] == time) )

    ax[3*columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['i_-'].where(data['Time'] == time) )

    ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
       data['i_Φ2'].where(data['Time'] == time) )

    ax[3*columns+2].plot(data['Position'].where(data['Time'] == time ),
             data['i_μe'].where(data['Time'] == time) )

    ax[2*columns+3].plot(data['Position'].where(data['Time'] == time ),
       data['dΦ_2/dx'].where(data['Time'] == time) )

    ax[3*columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['dμ_e/dx'].where(data['Time'] == time)/Fconst )


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[columns+1].set_title('Mass Average Velocity vs Time')
ax[columns+1].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('$\mathregular{N_+ \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+2].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+3].set_ylabel('$\mathregular{N_- \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2*columns+1].set_ylabel('$\mathregular{i_+ \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_ylabel('$\mathregular{i_- \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_xlabel('Position (cm)')

ax[2*columns+2].set_ylabel('$\mathregular{i_{Φ_2} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+2].set_ylabel('$\mathregular{i_{μ_e} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+2].set_xlabel('Position (cm)')

ax[2*columns+3].set_ylabel('$\mathregular{\\nabla Φ_2 \ (cm^{-1}) }$')
ax[3*columns+3].set_ylabel('$\mathregular{\\nabla μ_e \ (cm^{-1}) }$')
ax[3*columns+3].set_xlabel('Position (cm)')

fig.tight_layout()

# In[23]:
os.chdir(cw)
os.chdir('23_2_BinarySalt_ConcSoln_withCurrent')
# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

time_unique = np.unique(data['Time'])

print(data.keys())

rows, columns = [4, 3]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()


for time in [0, time_unique[1], 1, 10, 30, 300]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )

    ax[columns+2].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_+'].where(data['Time'] == time) )

    ax[columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_-'].where(data['Time'] == time) )


    ax[2*columns+1].plot(data['Position'].where(data['Time'] == time ),
       data['i_+'].where(data['Time'] == time) )

    ax[3*columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['i_-'].where(data['Time'] == time) )

    # ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
    #    data['i_Φ2'].where(data['Time'] == time) )
    # Phi_2 = data['Φ_2'].where(data['Time'] == time).dropna().values
    # print(Phi_2[-2], Phi_2[-1])
    ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
       data['Φ_2'].where(data['Time'] == time) )

    ax[3*columns+2].plot(data['Position'].where(data['Time'] == time ),
             data['i_μe'].where(data['Time'] == time) )

    ax[2*columns+3].plot(data['Position'].where(data['Time'] == time ),
       data['dΦ_2/dx'].where(data['Time'] == time) )

    ax[3*columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['dμ_e/dx'].where(data['Time'] == time)/Fconst )


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[columns+1].set_title('Mass Average Velocity vs Time')
ax[columns+1].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('$\mathregular{N_+ \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+2].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+3].set_ylabel('$\mathregular{N_- \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2*columns+1].set_ylabel('$\mathregular{i_+ \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_ylabel('$\mathregular{i_- \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_xlabel('Position (cm)')

# ax[2*columns+2].set_ylabel('$\mathregular{i_{Φ_2} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[2*columns+2].set_ylabel('$\mathregular{Φ_2}$')
ax[3*columns+2].set_ylabel('$\mathregular{i_{μ_e} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+2].set_xlabel('Position (cm)')

ax[2*columns+3].set_ylabel('$\mathregular{\\nabla Φ_2 \ (cm^{-1}) }$')
ax[3*columns+3].set_ylabel('$\mathregular{\\nabla μ_e \ (cm^{-1}) }$')
ax[3*columns+3].set_xlabel('Position (cm)')

fig.tight_layout()

# In[23]:
os.chdir(cw)
os.chdir('23_3_BinarySalt_ConcSoln_withCurrent')
# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

time_unique = np.unique(data['Time'])

print(data.keys())

rows, columns = [4, 3]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()


for time in [0, time_unique[1], 1, 10, 30, 300]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )

    ax[columns+2].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_+'].where(data['Time'] == time) )

    ax[columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_-'].where(data['Time'] == time) )


    ax[2*columns+1].plot(data['Position'].where(data['Time'] == time ),
       data['i_+'].where(data['Time'] == time) )

    ax[3*columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['i_-'].where(data['Time'] == time) )

    # ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
    #    data['i_Φ2'].where(data['Time'] == time) )
    ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
       data['Φ_2'].where(data['Time'] == time) )

    ax[3*columns+2].plot(data['Position'].where(data['Time'] == time ),
             data['i_μe'].where(data['Time'] == time) )

    ax[2*columns+3].plot(data['Position'].where(data['Time'] == time ),
       data['dΦ_2/dx'].where(data['Time'] == time) )

    ax[3*columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['dμ_e/dx'].where(data['Time'] == time)/Fconst )


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[columns+1].set_title('Mass Average Velocity vs Time')
ax[columns+1].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('$\mathregular{N_+ \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+2].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+3].set_ylabel('$\mathregular{N_- \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2*columns+1].set_ylabel('$\mathregular{i_+ \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_ylabel('$\mathregular{i_- \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_xlabel('Position (cm)')

# ax[2*columns+2].set_ylabel('$\mathregular{i_{Φ_2} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[2*columns+2].set_ylabel('$\mathregular{Φ_2}$')
ax[3*columns+2].set_ylabel('$\mathregular{i_{μ_e} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+2].set_xlabel('Position (cm)')

ax[2*columns+3].set_ylabel('$\mathregular{\\nabla Φ_2 \ (cm^{-1}) }$')
ax[3*columns+3].set_ylabel('$\mathregular{\\nabla μ_e \ (cm^{-1}) }$')
ax[3*columns+3].set_xlabel('Position (cm)')

fig.tight_layout()


# In[23]:
os.chdir(cw)
os.chdir('23_4_BinarySalt_ConcSoln_withCurrent')
# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

time_unique = np.unique(data['Time'])

print(data.keys())

rows, columns = [4, 3]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()


for time in [0, time_unique[1], 1, 10, 30, max(time_unique)]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )

    ax[columns+2].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_+'].where(data['Time'] == time) )

    ax[columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_-'].where(data['Time'] == time) )


    ax[2*columns+1].plot(data['Position'].where(data['Time'] == time ),
       data['i_+'].where(data['Time'] == time) )

    ax[3*columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['i_-'].where(data['Time'] == time) )

    # ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
    #    data['i_Φ2'].where(data['Time'] == time) )
    Phi_2 = data['Φ_2'].where(data['Time'] == time).dropna().values
    ax[2*columns+2].plot(data['Position'].where(data['Time'] == time ),
       data['Φ_2'].where(data['Time'] == time) - Phi_2[-1] )

    ax[3*columns+2].plot(data['Position'].where(data['Time'] == time ),
             data['i_μe'].where(data['Time'] == time) )

    ax[2*columns+3].plot(data['Position'].where(data['Time'] == time ),
       data['dΦ_2/dx'].where(data['Time'] == time) )

    ax[3*columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['dμ_e/dx'].where(data['Time'] == time)/Fconst )


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[columns+1].set_title('Mass Average Velocity vs Time')
ax[columns+1].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('$\mathregular{N_+ \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+2].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[columns+3].set_ylabel('$\mathregular{N_- \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2*columns+1].set_ylabel('$\mathregular{i_+ \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_ylabel('$\mathregular{i_- \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+1].set_xlabel('Position (cm)')

# ax[2*columns+2].set_ylabel('$\mathregular{i_{Φ_2} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[2*columns+2].set_ylabel('$\mathregular{Φ_2}$')
ax[3*columns+2].set_ylabel('$\mathregular{i_{μ_e} \ (A\ cm^{-2}\ s^{-1}) }$')
ax[3*columns+2].set_xlabel('Position (cm)')

ax[2*columns+3].set_ylabel('$\mathregular{\\nabla Φ_2 \ (cm^{-1}) }$')
ax[3*columns+3].set_ylabel('$\mathregular{\\nabla μ_e \ (cm^{-1}) }$')
ax[3*columns+3].set_xlabel('Position (cm)')

fig.tight_layout()


# In[24]:
os.chdir(cw)
os.chdir('24_BinarySalt_ConcSoln_withCurrent_Potential')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

time_unique = np.unique(data['Time'])
position = np.unique(data['Position'])

print(data.keys())

rows, columns = [2, 5]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()


for time in [0, time_unique[1], 1, 10, 30, 300]:
    ax[1].plot(data['Position'].where(data['Time'] == time ),
             data['ω_e'].where(data['Time'] == time ) )

    ax[columns+1].plot(data['Position'].where(data['Time'] == time ),
             data['velocity'].where(data['Time'] == time)*1e7 )

    ax[2].plot(data['Position'].where(data['Time'] == time ),
             data['flux_e'].where(data['Time'] == time) )

    ax[columns+2].plot(data['Position'].where(data['Time'] == time ),
              data['flux_0'].where(data['Time'] == time) )

    ax[3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_+'].where(data['Time'] == time) )

    ax[columns+3].plot(data['Position'].where(data['Time'] == time ),
             data['flux_-'].where(data['Time'] == time) )

    ax[4].plot(data['Position'].where(data['Time'] == time ),
          data['i_+'].where(data['Time'] == time) )

    ax[columns+4].plot(data['Position'].where(data['Time'] == time ),
             data['i_-'].where(data['Time'] == time) )

    ax[5].plot(data['Position'].where(data['Time'] == time ),
          data['Φ_2'].where(data['Time'] == time) )

    ax[columns+5].plot(data['Position'].where(data['Time'] == time ),
             data['dμ_e/dx'].where(data['Time'] == time)/Fconst )

    # ax[5].plot(data['Position'].where(data['Time'] == time ),
    #       data['dΦ_2/dx'].where(data['Time'] == time) )

    # dPhi2_ = data['dΦ_2/dx'].where(data['Time'] == time).dropna().values
    # Phi2 = np.zeros(len(dPhi2_))
    # for i in range(1, len(dPhi2_)):
    #     Phi2[i] = Phi2[i-1] + dPhi2_[i-1] * (position[i] - position[i-1])
    #
    # ax[5].plot(data['Position'].where(data['Time'] == time ).dropna().values,
    #       Phi2 )


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[columns+1].set_title('Mass Average Velocity vs Time')
ax[columns+1].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('$\mathregular{N_+ \ (mol\ cm^{-2}\ s^{-1}) }$')

ax[2].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[columns+1].set_xlabel('Position (cm)')

ax[columns+2].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[columns+2].set_xlabel('Position (cm)')

ax[columns+3].set_ylabel('$\mathregular{N_- \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[columns+3].set_xlabel('Position (cm)')

ax[4].set_ylabel('$\mathregular{i_+ \ (A\ cm^{-2}\ s^{-1}) }$')
ax[columns+4].set_ylabel('$\mathregular{i_- \ (A\ cm^{-2}\ s^{-1}) }$')
ax[columns+4].set_xlabel('Position (cm)')

fig.tight_layout()


# In[25]:
os.chdir(cw)
os.chdir('25_BinarySalt_ConcSoln_withCurrent_Potential')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [2, 3]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

z_cat = +1
z_an  = -1
nu_cat = 1.0
nu_an  = 1.0

data_slice = data.where(data['Position'] < max(data['Position']))
for time in [0, 1, 10, 30, 300]:
    ax[1].plot(data_slice['Position'].where(data['Time'] == time ),
             data_slice['ω_e'].where(data['Time'] == time ) )

    ax[2].plot(data_slice['Position'].where(data['Time'] == time ),
             data_slice['velocity'].where(data['Time'] == time)*1e7 )

    ax[3].plot(data_slice['Position'].where(data['Time'] == time ),
             data_slice['Φ_2'].where(data['Time'] == time) )

    N_e = data_slice['flux_e'].where(data['Time'] == time)
    N_cat = nu_cat * N_e
    N_an  = nu_an * N_e
    ax[4].plot(data_slice['Position'].where(data['Time'] == time ), N_e )

    N_0 = data_slice['flux_0'].where(data['Time'] == time)
    ax[5].plot(data_slice['Position'].where(data['Time'] == time ), N_0)


    i2 = Fconst * (z_cat * N_cat + z_an*N_an)
    print(z_cat * N_cat.dropna().values + z_an * N_an.dropna().values)
    ax[6].plot(data_slice['Position'].where(data['Time'] == time ), i2)


ax[1].set_title('Mass Fraction vs Time')
ax[1].set_ylabel('$\mathregular{ω_e}$')
# ax[1].set_xlabel('Position (cm)')

ax[2].set_title('Mass Average Velocity vs Time')
ax[2].set_ylabel('Velocity (nm/s)')

ax[3].set_ylabel('Φ_2')

ax[4].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[4].set_xlabel('Position (cm)')

ax[5].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')
ax[5].set_xlabel('Position (cm)')

ax[6].set_ylabel('i_2')
ax[6].set_xlabel('Position (cm)')

fig.tight_layout()

# In[1000]:
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

def Upwind(N,x0,xe,t0,te,m,c,alpha,f,boundary='Dirichlet'):

    # ----------------------
    ### Timestep

    t = np.linspace(t0, te, N)
    dt = t[1]-t[0]

    # ----------------------
    ### Direction conditions

    if c == 0:

        print('u(x,t) = f(x)')

    else:

        if c > 0:

            x = np.linspace(xe, x0, N)

        else:

            x = np.linspace(x0, xe, N)


        dx = x[1]-x[0]
        r = c*np.ones(N)*dt/dx

        # ----------------------
        ### Stability conditions

        if abs(r[0]) > 1+1e-15:

            print('r = {}'.format(r[0]))
            print('This solution is unstable, abs(r) must be less \
                  than or equal to one for stability.')

        # ----------------------
        ### Derivative matrix

        DL = spdiags([np.ones(N)], (0), N, N).tocsr() # NOT NECESSARY, MORE FOR
                                                      # CONTINUITY WITH OTHER SCHEMES
                                                      # THAT HAVE A LHS
        DR = spdiags([-r, 1+r], (-1, 0), N, N).tocsr()

        # ----------------------
        ### Boundary conditions

        if boundary == 'Dirichlet': # AGAIN, NOT NECESSARY, MORE FOR
                                    # CONTINUITY. HERE u(.,t) = 0.

            BL = spdiags([0], (0), N, N)
            BR = spdiags([0], (0), N, N)

        elif boundary == 'Neumann': # NEEDS TO BE FIXED

            BL = spdiags([0], (0), N, N)
            BR = spdiags([0], (0), N, N)

        else:

            BL = spdiags([0], (0), N, N)
            BR = spdiags([-r], (N-1), N, N)

        # ----------------------
        ### Initial condition

        u = f(x)

        # ----------------------
        ### Loop

        for k in range(m+1):

            RHS = (DR + BR)*u
            if boundary == 'Neumann':
                RHS[0] += c*dt*(alpha-u[0]/dx)
                # RHS[0] += c*dt*(alpha-(u[1] - u[0])/dx)

            u = spsolve(DL+BL, RHS)

            # u = spsolve(DL+BL, (DR+BR)*u)

            if k % 25 == 0:

                plt.plot(x,u)
# In[101]:
#       N,  x0,xe,t0,te,  m,c,alpha,f,                       boundary='Dirichlet'
Upwind(200,-40,40,0, 40, 100,-2,1,   lambda x: np.exp(-x**2), boundary='Periodic')

# In[101]:
Upwind(200,-40,40,0,40,100,-2,1,lambda x: np.exp(-x**2), boundary='Dirichlet')
# In[101]:
Upwind(200,-40,40,0,40,100,-2,0.1,lambda x: np.exp(-x**2), boundary='Neumann')

# In[25]:
os.chdir(cw)
os.chdir('30_UpWind')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
time_unique = np.unique(data['Time'])
print(data.keys())

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
# ax_t = ax[1].twinx()

z_cat = +1
z_an  = -1
nu_cat = 1.0
nu_an  = 1.0
Fconst = 96485

data_slice = data.where(data['Position'] < max(data['Position']))



for time in time_unique[::10]:#, 1, 10, 30, 300]:
    ax[1].plot(data_slice['Position'].where(data['Time'] == time ),
             data_slice['u'].where(data['Time'] == time ) )

#     ax[2].plot(data_slice['Position'].where(data['Time'] == time ),
#              data_slice['velocity'].where(data['Time'] == time)*1e7 )
#
#     ax[3].plot(data_slice['Position'].where(data['Time'] == time ),
#              data_slice['Φ_2'].where(data['Time'] == time) )
#
#     N_e = data_slice['flux_e'].where(data['Time'] == time)
#     N_cat = nu_cat * N_e
#     N_an  = nu_an * N_e
#     ax[4].plot(data_slice['Position'].where(data['Time'] == time ), N_e )
#
#     N_0 = data_slice['flux_0'].where(data['Time'] == time)
#     ax[5].plot(data_slice['Position'].where(data['Time'] == time ), N_0)
#
#
#     i2 = Fconst * (z_cat * N_cat + z_an*N_an)
#     print(z_cat * N_cat.dropna().values + z_an * N_an.dropna().values)
#     ax[6].plot(data_slice['Position'].where(data['Time'] == time ), i2)
#
#
# ax[1].set_title('Mass Fraction vs Time')
# ax[1].set_ylabel('$\mathregular{ω_e}$')
# # ax[1].set_xlabel('Position (cm)')
#
# ax[2].set_title('Mass Average Velocity vs Time')
# ax[2].set_ylabel('Velocity (nm/s)')
#
# ax[3].set_ylabel('Φ_2')
#
# ax[4].set_ylabel('$\mathregular{N_e \ (mol\ cm^{-2}\ s^{-1}) }$')
# ax[4].set_xlabel('Position (cm)')
#
# ax[5].set_ylabel('$\mathregular{N_0 \ (mol\ cm^{-2}\ s^{-1}) }$')
# ax[5].set_xlabel('Position (cm)')
#
# ax[6].set_ylabel('i_2')
# ax[6].set_xlabel('Position (cm)')

fig.tight_layout()
