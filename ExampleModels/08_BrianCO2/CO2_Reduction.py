'''
Plot Concentration Cell Data
LiPF6 in EC : EMC

Li-EMI-TFSI binary mixtures

'''

'''
Multi-Cursor
For a temporary fix:

    Open the console (Ctrl + Shift + I) (Alt + Cmd + I on MacOS)
    run atom.config.set('core.editor.multiCursorOnClick', true);

    Run all - hot-key (Atom - Hydrogen): Ctrl + Apple + Enter
'''

# In[0]:
# import libraries and helper functions
# Please see PythonHelperFunctions for what exactly is imported
import sys
sys.path.insert(0, '/Users/nicholasbrady/Documents/School/Academic/West Research/Projects/')
from PythonHelperFunctions import *
plot_parameters()
csfont = {'fontname':'Serif'}
# latfont = {'Computer Modern Roman'}
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# matplotlib.rc('text', usetex=False)
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# matplotlib.rc('text', usetex=True)
# https://matplotlib.org/users/customizing.html
matplotlib.rc('xtick', top=True, bottom=True, direction='in')
matplotlib.rc('ytick', left=True, right=True, direction='in')
plt.rc("axes.spines", top=True, right=True)
matplotlib.rc('axes', edgecolor='k')

import scipy
import timeit
from tqdm import tqdm
import math
# from sklearn.model_selection import train_test_split

def axes(number, rows, columns):
    fig = plt.figure(number, figsize=(6*columns, 4*rows), dpi=80)
    ax = []
    ax.append([])
    for i in range(1, rows*columns+1):
        ax.append(fig.add_subplot(rows, columns, i))

    return ax, fig

# os.getcwd()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/10_CO2_Reduction_SteadyState/')

# In[1]:
# In[20]:
from scipy.optimize import fsolve
'''
Solve the equilibrium system of equations for the paper Mathematical Modeling of CO2 Reduction to CO in Aqueous Electrolytes: I. Kinetic Study on Planar Silver and Gold Electrodes by Charles Delacourt, Paul L. Ridgway, and John Newman
'''
K_CO2 = 10**(-1.491)
P_CO2 = 1.013

K_pH = 10**(-13.74)
K_5  = 10**(6.05)
K_6  = 10**(9.73)

z_CO2  =  0.0
z_H    = +1.0
z_K    = +1.0
z_Na   = +1.0
z_OH   = -1.0
z_HCO3 = -1.0
z_CO3  = -2.0
z_ClO4 = -1.0

def equations(p):
    CO2, H, OH, HCO3, CO3, K = p

    return (CO2 - K_CO2 * P_CO2,
            K_6 * CO3 * H - HCO3,
            K_5 * H * HCO3 - CO2,
            K_pH - OH * H,
            K - 0.500,
            z_CO2 * CO2 + z_H * H + z_OH * OH + z_HCO3 * HCO3 + z_CO3 * CO3 + z_K * K)


CO2, H, OH, HCO3, CO3, K = fsolve(equations, (K_CO2 * P_CO2, 1e-7, 1e-7, 0.5, 0, 0.5))

print(CO2, H, OH, HCO3, CO3, K)
print(equations((CO2, H, OH, HCO3, CO3, K)))

# In[2]:
CER_equil = pd.read_table('Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
CER_equil.keys()


rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

CER_equil_slice = CER_equil[CER_equil['Position'] == 0.0]

ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[CO2]'], '-', color='blue')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[OH-]'], '-', color = 'darkorange')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[HCO3-]'], '-', color='green')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[CO3-2]'], '-', color = 'red')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[K+]'], '-', color = 'purple')
# ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[H+]'], '-')

# ax[1].set_xlabel('Electrode Potential, $\mathregular{\Phi_1}$ (V)')
ax[1].set_ylabel('Specie Concentration, [mol/L]')

ax[1].text(-0.55, 0.75, '$\mathregular{CO_2}$', color = 'blue', ha='right')
ax[1].text(-1.4, 0.31, '$\mathregular{OH^-}$', color = 'darkorange', ha='right')
ax[1].text(-1.25, 0.3, '$\mathregular{CO_3^{-2}}$', color = 'red', va='top')
ax[1].text(-0.9, 0.4, '$\mathregular{HCO_3^-}$', color = 'green', ha='right')
ax[1].text(-1.25, 0.67, '$\mathregular{K^+}$', color = 'purple')

ax[2].plot(CER_equil_slice['Phi_1'], -np.log10(CER_equil_slice['[H+]']), 'k-')

ax2_twin = ax[2].twinx()
ax2_twin.plot(CER_equil_slice['Phi_1'], CER_equil_slice['Phi_2']*1e3, 'r-')
ax2_twin.set_ylabel('Solution Potential, $\mathregular{\Phi_2}$ (mV)', color ='red')
ax2_twin.spines['right'].set_color('red')
ax2_twin.tick_params(axis='y', colors='red')
print(min(CER_equil_slice['Phi_2']*1e3))

ax[2].set_xlabel('Electrode Potential, $\mathregular{\Phi_1}$ (V)')
ax[2].set_ylabel('Solution pH')

fig.tight_layout()
fig.savefig('CER_Equilibrated.pdf', format='pdf', dpi=100, bbox_inches = "tight")

# In[3]:
rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

rows, columns = [1, 1]
ax2, fig2 = axes(2, rows, columns)

CER_equil_slice = CER_equil[CER_equil['Position'] == 0.0]

ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[CO2]'], '-', color='blue')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[OH-]'], '-', color='darkorange')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[HCO3-]'], '-', color='green')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[CO3-2]'], '-', color='red')
ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[K+]'], '-', color='purple')
# ax[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['[H+]'], '-')

ax[1].set_xlim([-2.0, -0.7])
ax[1].set_ylim([0.0, 0.8])
# ax[1].set_xlabel('Electrode Potential, $\mathregular{\Phi_1}$ (V)')
ax[1].set_ylabel('Specie Concentration, [mol/L]')

ax[1].text(-0.78, 0.2, '$\mathregular{CO_2}$', color = 'blue', ha='right')
ax[1].text(-1.42-0.103, 0.31, '$\mathregular{OH^-}$', color = 'darkorange', ha='right')
ax[1].text(-1.25-0.103, 0.3, '$\mathregular{CO_3^{-2}}$', color = 'red', va='top')
ax[1].text(-0.9-0.103, 0.4, '$\mathregular{HCO_3^-}$', color = 'green', ha='right')
ax[1].text(-1.25-0.103, 0.66, '$\mathregular{K^+}$', color = 'purple')

ax[2].plot(CER_equil_slice['Phi_1'], -np.log10(CER_equil_slice['[H+]']), '-')

ax[2].set_xlim([-2.0, -0.7])
ax[2].set_xlabel('Electrode Potential, $\mathregular{\Phi_1}$ (V)')
ax[2].set_ylabel('Solution pH')

ax2[1].plot(CER_equil_slice['Phi_1'], CER_equil_slice['RDE_Current']/1e3, '-', color = 'black')
ax2[1].set_xlabel('Electrode Potential, $\mathregular{\Phi_1}$ (V)')
ax2[1].set_ylabel('Current Density (A/$\mathregular{cm^2}$)')
ax2[1].set_xlim([-2.0, -0.7])
ax2[1].set_ylim([-0.15, 0.0])
ax2[1].tick_params(axis='both', which='major', pad=8)

ax2_twin = ax2[1].twinx()
di_dV = (CER_equil_slice['RDE_Current'].values[1:] - CER_equil_slice['RDE_Current'].values[:-1]) / (CER_equil_slice['Phi_1'].values[1:] - CER_equil_slice['Phi_1'].values[:-1])

ax2_twin.plot(CER_equil_slice['Phi_1'][:-1], di_dV/1e3, '-', color = 'red')
ax2_twin.spines['right'].set_color('red')
ax2_twin.tick_params(axis='y', colors='red')
ax2_twin.set_ylabel('di/dV (A/(V $\cdot \mathregular{cm^2}$))', color ='red')

fig.tight_layout()
