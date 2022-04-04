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
import scipy
from scipy import special

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

# print(os.getcwd())
cw = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/'
os.chdir(cw)

color_list = ['red', 'blue', 'green', 'darkorange', 'purple', 'black']

# In[1]:
# Diffusion Coefficient
rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

# conc range
c_Li = np.linspace(1e-5,1e-2, 200)
min_D_Li = 1e-7

# exponential decay
D_Li_0 = 1e-5
c_Li_0 = 9.0e-4
D_Li = 1.5*D_Li_0 * np.exp(-c_Li / c_Li_0) + min_D_Li
ax[1].semilogy(c_Li * 1e3, D_Li, 'k-')
ax[2].loglog(c_Li * 1e3, D_Li, 'k-')

# exponential subtraction
D_Li_0 = 1e-6
c_Li_0 = 5.0e-4
D_Li = 1e-5 - D_Li_0 * np.exp(c_Li*2.3e2)
ax[1].semilogy(c_Li * 1e3, D_Li, 'g-')
ax[2].loglog(c_Li * 1e3, D_Li, 'g-')

# erfc function
D_Li = scipy.special.erfc((c_Li - 0.01e-3)/20e-4)/2.0 * D_Li_0*10*2 + min_D_Li
ax[1].semilogy(c_Li * 1e3, D_Li, 'r-')
ax[2].loglog(c_Li * 1e3, D_Li, 'r-')

ax[1].set_xlabel("$c_{\mathregular{Li^+}}$ (mol / L)")
ax[2].set_xlabel("$c_{\mathregular{Li^+}}$ (mol / L)")

ax[1].set_ylabel("$D_{\mathregular{Li^+}}$ $\mathregular{(cm^2 / s)}$")
ax[2].set_ylabel("$D_{\mathregular{Li^+}}$ $\mathregular{(cm^2 / s)}$")

# In[2]:
os.chdir(cw)
os.chdir('11_BinarySalt_ILE_DiluteSoln')
rows, columns = [1, 3]
ax, fig = axes(1, rows, columns)

# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

for time in [0, 0.1, 1, 3, 10]:
    ax[1].plot(data['Position'].where(data['Time'] == time),
             data['Conc'].where(data['Time'] == time))
    # ax[2].plot(data['Position'].where(data['Time'] == time),
    #         data['Phi_1'].where(data['Time'] == time))

# plt.ylim(0, 1.1)
ax[1].set_title('Concentration Profiles vs Time')
ax[1].set_ylabel("$c_{\mathregular{Li^+}}$ (mol / L)")
ax[1].set_xlabel('Position (cm)')

df_NJ = data[data['Position'] == data['Position'].max()]
df_1 = data[data['Position'] == data['Position'].min()]

ax[2].plot(df_NJ['Time'], df_NJ['Phi_1']*1e3)
ax[2].plot(df_1['Time'], df_1['Phi_1']*1e3)

ax[2].set_ylabel('$\mathregular{\Phi_1}$ (mV)')
ax[2].set_xlabel('Time (hours)')

ax[3].plot(df_NJ['Time'], df_NJ['i_rxn']*1e3)

ax[3].set_ylabel('$\mathregular{i_{rxn} \ (mA/cm^2)}$ ')
ax[3].set_xlabel('Time (hours)')

fig.tight_layout()


# In[2]:
rows, columns = [2, 3]
ax, fig = axes(1, rows, columns)

for o, diff_coef in enumerate(['Constant/', 'Variable/']):
    os.chdir(cw)
    os.chdir('12_DiluteSoln_CyclicVolt_Steady/Rectangular/{}'.format(diff_coef))


    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    print(data.keys())

    ax[o*3+1].plot(data['DeltaV']*1e3, data['i_rxn']*1e3)
    ax[o*3+1].set_ylabel('$\mathregular{i_{rxn} \ (mA/cm^2)}$')

    df_NJ = data[data['Position'] == data['Position'].max()]
    df_1 = data[data['Position'] == data['Position'].min()]

    ax[o*3+2].plot(df_NJ['DeltaV']*1e3, df_NJ['Conc'])
    ax[o*3+2].plot(df_1['DeltaV']*1e3, df_1['Conc'])

    ax[o*3+2].set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')

    for delV in [min(abs(data['DeltaV'])), 50e-3, 100e-3, 150e-3, 400e-3]:
        data_ = data[data['DeltaV'] == delV]
        data_ = data_.sort_values('Position')
        ax[o*3+3].plot(data_['Position'], data_['Conc'])

    ax[o*3+3].set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')



ax[1].set_xticklabels([])
ax[2].set_xticklabels([])
ax[3].set_xticklabels([])

ax[4].set_ylabel('$\mathregular{i_{rxn} \ (mA/cm^2)}$')
ax[4].set_xlabel('$\mathregular{\Delta V}$ (mV)')

ax[5].set_xlabel('$\mathregular{\Delta V}$ (mV)')

ax[6].set_xlabel('Position (cm)')

fig.tight_layout()


# In[2]:
rows, columns = [2, 1]
fig = plt.figure(1, figsize=(10, 5*rows), dpi=80)

ax = []
ax.append([])
for i in range(1, rows*columns+1):
    ax.append(fig.add_subplot(rows, columns, i))




color_set = ['blue', 'orange', 'green', 'red', 'purple']

for o, diff_coef in enumerate(['Constant/', 'Variable/']):


    os.chdir(cw)
    os.chdir('12_DiluteSoln_CyclicVolt_Steady/Rectangular/{}'.format(diff_coef))


    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    print(data.keys())

    ax[o+1].plot(data['DeltaV']*1e3, data['i_rxn']*1e3, 'k-')

    # inset 1

    # if o == 0:
    #     axin1 = ax[o+1].inset_axes([0.08, 0.5, 0.35, 0.45])
    # elif o == 1:
    #     axin1 = ax[o+1].inset_axes([0.1, 0.5, 0.35, 0.45])

    axin1 = ax[o+1].inset_axes([0.08, 0.48, 0.35, 0.48])
    axin2 = ax[o+1].inset_axes([0.62, 0.15, 0.35, 0.48])

    df_NJ = data[data['Position'] == data['Position'].max()]
    df_1 = data[data['Position'] == data['Position'].min()]

    for oo, delV in enumerate([min(abs(data['DeltaV'])), 50e-3, 100e-3, 150e-3, 400e-3]):
        color = color_set[oo]

        data_ = data[data['DeltaV'] == delV]
        data_ = data_.sort_values('Position')
        axin2.plot(data_['Position'], data_['Conc'], color = color)
        ax[o+1].plot(delV*1e3, data_['i_rxn'].values[0]*1e3, 'o', color = color, markersize=10, markeredgecolor='k')

        axin2.text(0.5, 0.8, '$\mathregular{i_{rxn} > 0}$', fontsize = 20, transform = axin2.transAxes)


        data_ = data[data['DeltaV'] == -delV]
        data_ = data_.sort_values('Position')
        axin1.plot(data_['Position'], data_['Conc'], color = color)
        ax[o+1].plot(-delV*1e3, data_['i_rxn'].values[0]*1e3, 'o', color = color, markersize=10, markeredgecolor='k')

        axin1.text(0.1, 0.8, '$\mathregular{i_{rxn} < 0}$', fontsize = 20, transform = axin1.transAxes)
#
    axin1.set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')
    axin1.set_xlabel('Position (cm)')

    axin2.set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')
    axin2.set_xlabel('Position (cm)')

    axin1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.25))
    axin2.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.25))
    if o == 0:
        axin1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        axin1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.5))
        axin2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        axin2.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.5))

    ax[o+1].set_ylabel('$\mathregular{i_{rxn} \ (mA/cm^2)}$')

ax[1].set_xticklabels([])
ax[2].set_xlabel('$\mathregular{\Delta V}$ (mV)')

fig.tight_layout()

# In[2]:
rows, columns = [2, 1]
fig = plt.figure(1, figsize=(10, 5*rows), dpi=80)

ax = []
ax.append([])
for i in range(1, rows*columns+1):
    ax.append(fig.add_subplot(rows, columns, i))




color_set = ['blue', 'orange', 'green', 'red', 'purple']

for o, diff_coef in enumerate(['Constant/', 'Variable/']):


    os.chdir(cw)
    os.chdir('12_DiluteSoln_CyclicVolt_Steady/Spherical/{}'.format(diff_coef))


    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    print(data.keys())

    ax[o+1].plot(data['DeltaV']*1e3, data['i_rxn']*1e3, 'k-')

    # inset 1

    if o == 0:
        axin1 = ax[o+1].inset_axes([0.08, 0.5, 0.35, 0.45])
    elif o == 1:
        axin1 = ax[o+1].inset_axes([0.1, 0.5, 0.35, 0.45])

    axin2 = ax[o+1].inset_axes([0.62, 0.18, 0.35, 0.45])

    df_NJ = data[data['Position'] == data['Position'].max()]
    df_1 = data[data['Position'] == data['Position'].min()]

    for oo, delV in enumerate([min(abs(data['DeltaV'])), 50e-3, 100e-3, 150e-3, 400e-3]):
        color = color_set[oo]

        data_ = data[data['DeltaV'] == delV]
        data_ = data_.sort_values('Position')
        axin2.plot(data_['Position'], data_['Conc'], color = color)
        ax[o+1].plot(delV*1e3, data_['i_rxn'].values[0]*1e3, 'o', color = color, markersize=10, markeredgecolor='k')

        axin2.text(0.5, 0.8, '$\mathregular{i_{rxn} > 0}$', fontsize = 20, transform = axin2.transAxes)

        data_ = data[data['DeltaV'] == -delV]
        data_ = data_.sort_values('Position')
        axin1.plot(data_['Position'], data_['Conc'], color = color)
        ax[o+1].plot(-delV*1e3, data_['i_rxn'].values[0]*1e3, 'o', color = color, markersize=10, markeredgecolor='k')

        axin1.text(0.5, 0.1, '$\mathregular{i_{rxn} < 0}$', fontsize = 20, transform = axin1.transAxes)
#
    axin1.set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')
    axin1.set_xlabel('Position (cm)')

    axin2.set_ylabel('$c_{\mathregular{Li^+}}$ (mol/L)')
    axin2.set_xlabel('Position (cm)')


    ax[o+1].set_ylabel('$\mathregular{i_{rxn} \ (mA/cm^2)}$')

ax[1].set_xticklabels([])
ax[2].set_xlabel('$\mathregular{\Delta V}$ (mV)')

fig.tight_layout()


# In[4]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']

os.chdir(cw)
os.chdir('12_DiluteSoln_CyclicVolt_Steady/Rectangular/Variable/')


# import the data file
data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

data_ = data[data['DeltaV'] == 70e-3]

ax[1].plot(data_['Position'], data_['Conc'], 'k-')

data_ = data[data['DeltaV'] == 72e-3]

ax[1].plot(data_['Position'], data_['Conc'], 'k--')

# In[4]:
rows, columns = [1, 1]
fig = plt.figure(1, figsize=(10, 5*rows), dpi=80)

ax = []
ax.append([])
for i in range(1, rows*columns+1):
    ax.append(fig.add_subplot(rows, columns, i))

linestyle = ['-', '--']
lineColor = ['red', 'blue']

for oo, geo in enumerate(['Rectangular/', 'Spherical/']):

    for o, diff_coef in enumerate(['Constant/', 'Variable/']):
        os.chdir(cw)
        os.chdir('12_DiluteSoln_CyclicVolt_Steady/{}{}'.format(geo, diff_coef))


        ls = linestyle[o]
        lC = lineColor[oo]

        # import the data file
        data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

        data = data.sort_values('DeltaV')

        ax[1].plot(data['DeltaV']*1e3, abs(data['i_rxn'])*1e3, color = lC, linestyle=ls)

ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

# In[4]:
rows, columns = [1, 1]
fig = plt.figure(1, figsize=(10, 5*rows), dpi=80)

ax = []
ax.append([])
for i in range(1, rows*columns+1):
    ax.append(fig.add_subplot(rows, columns, i))

linestyle = ['-', '--']
lineColor = ['red', 'blue']

for oo, geo in enumerate(['Rectangular/', 'Spherical/']):

    for o, diff_coef in enumerate(['Constant/', 'Variable/']):
        os.chdir(cw)
        os.chdir('13_DiluteSoln_CyclicVolt_Steady_SatConc/{}{}'.format(geo, diff_coef))


        ls = linestyle[o]
        lC = lineColor[oo]

        # import the data file
        data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

        data = data.sort_values('DeltaV')

        ax[1].semilogy(data['DeltaV']*1e3, abs(data['i_rxn'])*1e3, color = lC, linestyle=ls)

ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

# In[10]:

'''
    Half-Cell Simulations
'''
os.chdir(cw)
HalfCell_dir = '15_DiluteSoln_CyclicVolt_HalfCell/'
os.chdir(HalfCell_dir)
HalfCell_dir = os.getcwd()


rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']
lineColor = ['blue', 'red']



for o, diff_coef in enumerate(['Constant/', 'Variable/']):
    os.chdir(HalfCell_dir)
    os.chdir(diff_coef)


    ls = '-'
    lC = lineColor[o]

    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    data_electrode = data[data['Position'] == data['Position'].min()]

    ax[1].plot(data['DeltaV']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    ax[1].plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)
    # ax[2].semilogy(data_electrode['DeltaV']*1e3, data_electrode['Conc'], color = lC, linestyle=ls)


ax[1].text(70, 5, 'Kinetic \nLimitation', ha='center', color='black')
ax[1].plot([90, 127], [4.9, 4.3], 'k-', linewidth = 1)

# i_lim_MT
Fconst = 96485
i_lim_MT = -1e-5 * 1e-3/1.0 * Fconst * 1e3

ax[1].plot(data['DeltaV']*1e3, np.ones(len(data))*i_lim_MT, 'b--', linewidth=2)
ax[1].set_ylim([-2, 6])
ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")



ax_inset = ax[1].inset_axes([0.19, 0.45, 0.5, 0.5])
# conc range
c_Li = np.logspace(-5,np.log10(2e-2), 200)
min_D_Li = 1e-7
D_Li_0 = 1e-5
c_Li_0 = 9.0e-4
D_Li = 1.5*D_Li_0 * np.exp(-c_Li / c_Li_0) + min_D_Li
ax_inset.loglog(c_Li * 1e3, D_Li, 'r-')
ax_inset.loglog(c_Li * 1e3, np.ones(len(c_Li))*D_Li_0, 'b-')

ax_inset.set_xlabel("$c_{\mathregular{Li^+}}$ (mol / L)", labelpad=-1)
ax_inset.set_ylabel("$D_{\mathregular{Li^+}}$ $\mathregular{(cm^2 / s)}$", labelpad=-3)

fig.savefig(HalfCell_dir + '/Variable_D0.pdf',
            format='pdf', dpi=300, bbox_inches = "tight")




# In[20]:
'''
    Half-Cell Simulations
'''
os.chdir(cw)
HalfCell_dir = '15_DiluteSoln_CyclicVolt_HalfCell/'
os.chdir(HalfCell_dir)
HalfCell_dir = os.getcwd()
os.chdir('Constant')
os.chdir('Vary_Diff0/')
cwd_vary_diff = os.getcwd()

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']
lineColor = ['blue', 'green']

dirs = [name for name in os.listdir('.') if os.path.isdir(name)]

for o, d0 in enumerate(dirs):
    os.chdir(cwd_vary_diff)
    os.chdir(d0)

    ls = '-'
    lC = color_list[o]

    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    data_electrode = data[data['Position'] == data['Position'].min()]

    ax[1].plot(data['DeltaV']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    ax[1].plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)

    Fconst = 96485
    i_lim_MT = -float(d0) * 1e-3/1.0 * Fconst * 1e3
    print(i_lim_MT)
    ax[1].plot(data['DeltaV']*1e3, np.ones(len(data))*i_lim_MT, '--', color = lC, linewidth = 2)

    # data_ = data[data['DeltaV'] == 100e-3]
    # data_ = data_.sort_values('Position')
    # ax[2].plot(data_['Position'], data_['Conc'], color = lC)
    #
    # data_ = data[data['DeltaV'] == -300e-3]
    # data_ = data_.sort_values('Position')
    # ax[3].plot(data_['Position'], data_['Conc'], color = lC)
    # ax[2].semilogy(data_electrode['DeltaV']*1e3, data_electrode['Conc'], color = lC, linestyle=ls)
#
#

ax[1].text(50, 4.5, 'Kinetic \nLimitation', ha='center', color='black')
ax[1].plot([90, 120], [4.4, 3.9], 'k-', linewidth = 1)

ax[1].text(-300, 0.7, 'Mass Transfer \nLimitations', ha='center', color='black')
ax[1].plot([-399, -360], [-1.75, 0.5], 'b-', linewidth = 1.5)
ax[1].plot([-350, -330], [-0.8, 0.5], 'r-', linewidth = 1.5)
ax[1].plot([-300, -300], [0.05, 0.5], 'g-', linewidth = 1.5)

ax[1].text(-270, -0.05, '$D_0 = 1 \cdot 10^{-6}$', color = 'green', ha='left', va = 'bottom', fontsize=12)
ax[1].text(-330, -0.9,  '$D_0 = 1 \cdot 10^{-5}$', color = 'red',  ha='left', va = 'bottom', fontsize=12)
ax[1].text(-380, -1.85,  '$D_0 = 2 \cdot 10^{-5}$', color = 'blue', ha='left', va = 'bottom', fontsize=12)

ax[1].set_xlim([-420, 199])
ax[1].set_ylim([-3, 6])
ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

fig.savefig(HalfCell_dir + '/Vary_Diff0.pdf',
            format='pdf', dpi=300, bbox_inches = "tight")



# In[20]:
'''
    Half-Cell Simulations
'''
os.chdir(cw)
HalfCell_dir = '15_DiluteSoln_CyclicVolt_HalfCell/'
os.chdir(HalfCell_dir)
HalfCell_dir = os.getcwd()
os.chdir('Constant')
os.chdir('Vary_Diff0/')
cwd_vary_diff = os.getcwd()

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']
lineColor = ['blue', 'green']

dirs = [name for name in os.listdir('.') if os.path.isdir(name)]
print(dirs)

lineColor

# ax_inset = ax[1].inset_axes([0.18, 0.62, 0.5, 0.33])

for o, d0 in enumerate(dirs):
    os.chdir(cwd_vary_diff)
    os.chdir(d0)

    ls = '-'
    lC = color_list[o]

    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    data_electrode = data[data['Position'] == data['Position'].min()]

    ax[1].plot(data['DeltaV']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    ax[1].plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)

    Fconst = 96485
    i_lim_MT = -float(d0) * 1e-3/1.0 * Fconst * 1e3
    print(i_lim_MT)
    ax[1].plot(data['DeltaV']*1e3, np.ones(len(data))*i_lim_MT, '--', color = lC, linewidth = 2)
    # ax[1].plot(data['DeltaV']*1e3, (data['i_kappa'])*1e3, color = 'darkorange', linestyle='--')

    # ax_inset.plot(data['Phi_2']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    # ax_inset.plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)

# ax_inset.plot(data['DeltaV']*1e3, (data['i_kappa'])*1e3, color = 'darkorange', linestyle='--')
# ax_inset.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(50))
# ax_inset.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
# # ax_inset.set_xlim([0, 400])
# # ax_inset.set_ylim([0, 21])
# ax_inset.set_xlabel("$\mathregular{\Delta V}$ (mV)", fontsize=13)
# ax_inset.set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$", labelpad=0, fontsize=13)

ax[1].text(50, 4.5, 'Kinetic \nLimitation', ha='center', color='black')
ax[1].plot([90, 120], [4.4, 3.9], 'k-', linewidth = 1)

ax[1].text(-390, 0.7, 'Mass Transfer Limitations', ha='left', color='black')
ax[1].plot([-399, -360], [-1.75, 0.5], 'b-', linewidth = 1.5)
ax[1].plot([-350, -330], [-0.8, 0.5], 'r-', linewidth = 1.5)
ax[1].plot([-300, -300], [0.05, 0.5], 'g-', linewidth = 1.5)

ax[1].text(-270, -0.05, '$D_0 = 1 \cdot 10^{-6}$', color = 'green', ha='left', va = 'bottom', fontsize=12)
ax[1].text(-330, -0.9,  '$D_0 = 1 \cdot 10^{-5}$', color = 'red',   ha='left', va = 'bottom', fontsize=12)
ax[1].text(-380, -1.9,  '$D_0 = 2 \cdot 10^{-5}$', color = 'blue',  ha='left', va = 'bottom', fontsize=12)

ax[1].set_xlim([-420, 150])
ax[1].set_ylim([-3, 6])
ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

fig.savefig(HalfCell_dir + '/Vary_Diff0.pdf',
            format='pdf', dpi=300, bbox_inches = "tight")

# In[20]:
'''
    Half-Cell Simulations
'''
os.chdir(HalfCell_dir)
os.chdir('Constant')
os.chdir('C_sat/')
cwd_vary_diff = os.getcwd()

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']
lineColor = ['blue', 'green']

dirs = [name for name in os.listdir('.') if os.path.isdir(name)]

for o, d0 in enumerate(dirs):
    os.chdir(cwd_vary_diff)
    os.chdir(d0)

    ls = '-'
    lC = color_list[o]

    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    data_electrode = data[data['Position'] == data['Position'].min()]

    ax[1].plot(data['DeltaV']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    ax[1].plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)

    # data_ = data[data['DeltaV'] == 100e-3]
    # data_ = data_.sort_values('Position')
    # ax[2].plot(data_['Position'], data_['Conc'], color = lC)
    #
    # data_ = data[data['DeltaV'] == -300e-3]
    # data_ = data_.sort_values('Position')
    # ax[3].plot(data_['Position'], data_['Conc'], color = lC)
    # ax[2].semilogy(data_electrode['DeltaV']*1e3, data_electrode['Conc'], color = lC, linestyle=ls)
#
#
# # i_lim_MT
# Fconst = 96485
# i_lim_MT = -1e-5 * 1e-3/1.0 * Fconst * 1e3
# print(i_lim_MT)
# ax[1].plot(data['DeltaV']*1e3, np.ones(len(data))*i_lim_MT, 'k--')

ax[1].text(-50, 1.0, 'Kinetic \nLimitation', ha='center', color='black')
ax[1].plot([40, -30], [0.7, 0.95], 'k-', linewidth = 1)

ax[1].text(120, 1.0,  '$c_{sat} = 2.0$',   color = 'blue',  ha='left', va = 'bottom', fontsize=17)
ax[1].text(120, 0.5,  '$c_{sat} = 1.5$', color = 'green', ha='left', va = 'bottom', fontsize=17)
ax[1].text(120, 0.15, '$c_{sat} = 1.2$', color = 'red',   ha='left', va = 'top',    fontsize=17)

ax[1].set_ylim([-1.2, 1.5])
ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

fig.savefig(HalfCell_dir + '/Vary_Csat.pdf',
            format='pdf', dpi=300, bbox_inches = "tight")

# In[20]:
'''
    Half-Cell Simulations
'''
os.chdir(HalfCell_dir)
os.chdir('Constant')
os.chdir('Vary_kappa/')
cwd_vary_diff = os.getcwd()

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

linestyle = ['-', '--']
lineColor = ['blue', 'green']

dirs = [name for name in os.listdir('.') if os.path.isdir(name)]

for o, d0 in enumerate(dirs):
    os.chdir(cwd_vary_diff)
    os.chdir(d0)

    ls = '-'
    lC = color_list[o]

    # import the data file
    data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

    data = data.sort_values('DeltaV')

    data_electrode = data[data['Position'] == data['Position'].min()]

    ax[1].plot(data['DeltaV']*1e3, (data['i_rxn'])*1e3, color = lC, linestyle=ls)
    ax[1].plot(data['DeltaV']*1e3, (data['i_BV'])*1e3,  color = 'black', linestyle='--', linewidth=3)

    data_ = data[data['DeltaV']*1e3 >= 800]
    p = np.polyfit(data_['DeltaV']*1e3, data_['i_rxn']*1e3, 1)
    delV = np.linspace(300, 1000)
    ax[1].plot(delV, np.polyval(p, delV), '--', color = 'k', linewidth=1.5)
    print(p)

ax[1].text(1000, 28,  '$\kappa = 1 \cdot 10^{-1} \ \mathregular{S/cm}$',
    color = 'green', ha='right', va = 'bottom', fontsize=17, rotation = 0)
ax[1].text(1000, 15,  '$\kappa = 3 \cdot 10^{-2}\  \mathregular{S/cm}$',
    color = 'blue',  ha='right', va = 'center', fontsize=17, rotation = 30)
ax[1].text(1000, 4, '$\kappa = 1 \cdot 10^{-2} \ \mathregular{S/cm}$',
    color = 'red',   ha='right', va = 'center',    fontsize=17, rotation = 10)

ax[1].set_ylim([0, 40])
ax[1].set_xlabel("$\mathregular{\Delta V}$ (mV)")
ax[1].set_ylabel("$\mathregular{i_{rxn}} \ \mathregular{(mA/cm^2)}$ ")

fig.savefig(HalfCell_dir + '/Vary_Kappa.pdf',
            format='pdf', dpi=300, bbox_inches = "tight")
