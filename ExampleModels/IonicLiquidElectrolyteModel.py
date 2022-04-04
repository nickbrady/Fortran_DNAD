'''
Written by Nicholas Brady
March 16, 2021

Run all - hot-key (Atom - Hydrogen): Ctrl + Apple + Enter
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

Rigc   = 8.314             # Ideal gas constant [J/(mol*K)]
Temp   = 298               # Temperature [K]
Fconst = 96485             # Faraday's Constant [C/mol]

z_1_Li    =  1.0           # cation charge Li^+
z_2_EMI   =  1.0           # cation charge, EMI^+
z_3_TFSI  = -1.0           # anion charge of TFSI^-

nu_A_Li   = 1.0
nu_A_TFSI = 1.0
nu_A      = nu_A_Li + nu_A_TFSI
nu_B_EMI  = 1.0
nu_B_TFSI = 1.0
nu_B      = nu_B_EMI + nu_B_TFSI

MW_A_LiTFSI  = 287.09     # molecular weight Li-TFSI  g/mol
MW_B_EMITFSI = 391.31     # molecular weight EMI-TFSI g/mol
MW_Li = 6.941
MW_TFSI = MW_A_LiTFSI - MW_Li * nu_A_Li
MW_EMI = MW_B_EMITFSI - MW_TFSI * nu_B_TFSI

diff_12_Li_EMI = 1e-7                       # diffusion coefficient [cm^2/s]
diff_13_Li_TFSI = 2e-7
diff_23_EMI_TFSI = 3e-7

applied_current_A = -1.0e-4 *5

def density(cA):
    return 1.65 - (1.5 - 1.65)/2e-3 * (cA - 2e-3)

def cA_to_xA(c_A):
    rho = density(c_A)
    c_B = (rho - MW_A_LiTFSI * c_A) / MW_B_EMITFSI
    x_A = c_A / (c_A + c_B)

    return x_A

def chemical_potential(xA):
    return Rigc * Temp * log(xA ** nu_A_Li)

def fundamental_diff(cA):
    rho = density(cA)    # rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB
    c_Total = c1 + c2 + c3

    diff_fund = z_3_TFSI**2 * c_Total / (nu_A_Li * nu_B_EMI) / \
                (z_3_TFSI**2 * c3 / diff_12_Li_EMI
                + z_2_EMI**2 * c2 / diff_13_Li_TFSI
                +  z_1_Li**2 * c1 / diff_23_EMI_TFSI)

    return diff_fund

def practical_diff(cA):
    rho = density(cA)    # rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI
    c3  = nu_A_TFSI * cA + nu_B_TFSI * cB

    # need to add activity relationship
    diff_prac = fundamental_diff(cA) * (c3 / (cA + cB)) * (nu_A_Li)
    return diff_prac

def transference_number_common_ion(cA):
    rho = density(cA)    # rho = MW_A * cA + MW_B * cB
    cB  = (rho - MW_A_LiTFSI * cA) / MW_B_EMITFSI

    c1      = nu_A_Li * cA
    c2      = nu_B_EMI * cB
    c3      = nu_A_TFSI * cA + nu_B_TFSI * cB

    t_1_c = z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23_EMI_TFSI - z_3_TFSI /diff_12_Li_EMI) / \
    ( (z_2_EMI / diff_13_Li_TFSI - z_3_TFSI / diff_12_Li_EMI) + z_1_Li * c1/ (z_2_EMI * c2) * (z_1_Li / diff_23_EMI_TFSI - z_3_TFSI / diff_12_Li_EMI) )

    return t_1_c





cA = 0.001
print(fundamental_diff(cA))
print(transference_number_common_ion(cA))

# In[2]:
os.chdir(cw)
os.chdir('7_BinarySaltILE')

# import the data file
data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

print(data.keys())

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

NJ = len(np.unique(data['Position']))

times_array = np.array([0.0, 0.1, 0.2, 0.5, 0.99, 1.0, 4.9, 5.0, 5.1, 5.2, 5.5, 6.0, 9.9])

for time in np.unique(data['Time']):
    if any(np.isclose(times_array - time, 0.0, atol=5e-3)):

        if time < 5.0:

            df_slice = data[(data['State'] == 'D') & (data['Time'] == time)]
            ax[1].plot(df_slice['Position'], df_slice['Soln_Conc'])

        if time >= 5.0:
            df_slice = data[(data['State'] == 'R') & (data['Time'] == time)]
            ax[2].plot(df_slice['Position'], df_slice['Soln_Conc'])


ax[1].set_ylabel('$c_A$  (mol/L)')
ax[1].set_xlabel('Position (μm)')

ax[2].set_ylabel('$c_A$  (mol/L)')
ax[2].set_xlabel('Position (μm)')

df_slice = data[(data['State'] == 'D') & (data['Position'] == data['Position'].max())]
ax[3].plot((df_slice['Time'])**0.5, df_slice['Soln_Conc'])
ax[3].set_xlim(0, 1)




end_discharge = df_slice['Time'].max()
df_slice = data[(data['State'] == 'R') & (data['Position'] == data['Position'].max())]
ax[4].plot(df_slice['Time'] - end_discharge, np.log(df_slice['Soln_Conc'] - 1.0))

df_slice_fit = df_slice[(df_slice['Time'] >= (0.5 + end_discharge)) &
                        (df_slice['Time'] <= (5.5 + end_discharge))]

ax[4].set_xlim([0, 6])

# fit slope of relaxation concentration
P = np.polyfit(df_slice_fit['Time'] - end_discharge, np.log(df_slice_fit['Soln_Conc'] - 1.0), 1)
relax_time = df_slice_fit['Time'] - end_discharge

ax[4].plot(relax_time, np.polyval(P, relax_time), 'r--', linewidth=3)
diff = -P[0]/3600. * (data['Position'].max() * 1e-4)**2 / np.pi**2
print(P[0], diff)


fig.tight_layout()

# In[5]:
data['Time']

# In[2]:
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)

times_array = np.linspace(0,1./60,5)
for time in np.unique(data['Time']):
    if any(np.isclose(times_array - time, 0.0, atol=5e-5)):
        df_slice = data[(data['State'] == 'D') & (data['Time'] == time)]
        ax[1].plot(df_slice['Position'], df_slice['Soln_Conc'])


ax[1].set_ylabel('$c_A$  (mol/L)')
ax[1].set_xlabel('Position (μm)')
ax[1].set_xlim(4800,5000)

df_slice = data[(data['State'] == 'D') & (data['Position'] == data['Position'].max())]
ax[2].plot((df_slice['Time']*3600)**0.5, df_slice['Soln_Conc'])
# ax[2].set_xlim(0, 1.)
ax[2].set_ylim(1.0, 1.12)
ax[2].set_xlabel('$\mathregular{t^{1/2}}$ ($\mathregular{s^{1/2}}$)')

cA = 0.001
D  = practical_diff(cA)
t_1_c = transference_number_common_ion(cA)

print(cA*1e3, D, t_1_c)

c = (1. - t_1_c) / (z_1_Li * nu_A_Li * Fconst) * np.sqrt(np.pi/D) * -applied_current_A * (df_slice['Time']*3600)**0.5 + cA
c = c.values * 1e3
ax[2].plot((df_slice['Time']*3600)**0.5, c)

c2 = 2*(1. - t_1_c) / (z_1_Li * nu_A_Li * Fconst) * np.sqrt(1/(np.pi*D)) * -applied_current_A * (df_slice['Time']*3600)**0.5 + cA
c2 = c2.values * 1e3
ax[2].plot((df_slice['Time']*3600)**0.5, c2)

# end_exp_c = df_slice['Soln_Conc'].values[-1] / 1e3
# D = practical_diff(end_exp_c)
# t = transference_number_common_ion(end_exp_c)
#
# print(end_exp_c*1e3, D, t)

# c3 = 2*(1. - t_1_c) / (z_1_Li * nu_A_Li * Fconst) * np.sqrt(1/(np.pi*D)) * -applied_current_A * (df_slice['Time']*3600)**0.5 + cA
# c3 = c3.values * 1e3
# ax[2].plot((df_slice['Time']*3600)**0.5, c3)


# In[5]:
import scipy.special

t = 60 # seconds
print(60/3600)

x = np.linspace(0, 200e-4, 1000)
eta = x/(2 * np.sqrt(D * t))

c = cA + applied_current_A * ( 1. - t_1_c) / (z_1_Li * nu_A_Li * Fconst * np.sqrt(D)) * (2 * np.sqrt(t/np.pi) * np.exp(-x**2/(4*D*t)) - x/np.sqrt(D) * scipy.special.erfc(x/(2*np.sqrt(D*t))) )

plt.plot(x * 1e4, c * 1e3)

df_slice = data[(data['State'] == 'D') & (data['Time'] == 0.016640000000000002)]
plt.plot(df_slice['Position'], df_slice['Soln_Conc'])
plt.xlim(0, 200)
plt.ylim(0.85, 1.01)
