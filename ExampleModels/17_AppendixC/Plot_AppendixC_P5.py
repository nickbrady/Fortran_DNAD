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
from sklearn.model_selection import train_test_split


def axes(number, rows, columns):
    fig = plt.figure(number, figsize=(6*columns, 5*rows), dpi=80)
    ax = []
    ax.append([])
    for i in range(1, rows*columns+1):
        ax.append(fig.add_subplot(rows, columns, i))

    return ax, fig

colors = ['blue', 'darkorange', 'green', 'purple']
markers = ['o', 's', '^']

Rigc            = 8.314
Fconst          = 96485
Temp            = 298
alpha_a         = 0.5
alpha_c         = 0.5

const = 3
a =  0.51023262
b = -0.61592201


D = 1e-6        # cm2/s
nu = 1e-2       # cm2/s

# In[00]:
'''
    Rotating Disk Electrode (RDE)
    Chapter 15.4, 17.2, 19.1
'''


# In[1]:
'''
Velocity Profile Approximations

Case large Schmidt number,  Sc > 100
small values of ζ

Equation 15.29

H = -aζ^2 + 1/3*ζ^3 + b/6*ζ^4 + ...
'''

zeta = np.linspace(0,1.2)

H = -a * zeta**2 + 1/3 * zeta**3 + b/6 * zeta**4

plt.plot(zeta, H, 'blue')

'''
Equation 19.2

H = -aζ^2
'''

plt.plot(zeta, -a *zeta**2, 'k--')

'''
    Equation from Appendix C - PROGRAM MIGR (extension of equation 15.29)

    H = -aζ^2 + 1/3*ζ^3 + b/6*ζ^4 + b^2/30*ζ^5 + 3a/180*ζ^6 - (1 - 4ab)/1260*ζ^7 + ...
'''

zeta = np.linspace(0,2)

H = -a*zeta**2*(1. - zeta/a*(1/3 + b/6*zeta + (b+zeta)**2/30 + a*zeta*3/180 - (1 - 4*a*b)*zeta**4/1260))
#
plt.plot(zeta, H , 'red')


'''
Case small Schmidt number, Sc < 0.1
large values of ζ

Equation 15.31

H = -α + 2A/α * exp(-α*ζ) + ...
'''

zeta = np.linspace(1,4)
alpha = 0.88447
A = 0.934

H = -alpha + 2*A/alpha * np.exp(-alpha*(zeta)) # +0.35

plt.plot(zeta, H, 'darkorange')

plt.xlabel('ζ = $z \sqrt{Ω/ν}$')
plt.ylabel('H')

# In[2]:
'''
    Chapter 15.4 - reproduce figure 15.2
    Von Karman Approximation for a rotating disk (infinite size)

    Solve equations                 | x = 0                 | x = ∞
    2F + H' = 0                     | F = 0                 | F = 0
    F^2 - G^2 + HF' = F''           | H = 0                 | 2F + H' = 0  (--> H' = 0)
    2FG + HG' = G''                 | G = 1                 | G = 0
'''
ch15_4_dir = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/Ch15_4'
os.chdir(ch15_4_dir)

velocity_data = pd.read_table('Velocity_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

it = velocity_data['Time'].max()

velocity_data = velocity_data[velocity_data['Time'] == it]


rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

ax[1].plot(velocity_data['F'], velocity_data['Position'], color='blue')
ax[1].plot(velocity_data['G'], velocity_data['Position'], color='darkorange')
ax[1].plot(-velocity_data['H'], velocity_data['Position'], color='green')

ax[1].set_ylim([0, 4])
ax[1].set_xlim([0, 1])

ax[1].text(0.15, 1.25, 'F', ha='right', color='blue')
ax[1].text(0.6, 0.8, 'G', ha='left', color='darkorange')
ax[1].text(0.73, 3, '$-$H', ha='right', color='green')

ax[1].set_ylabel('$\mathregular{ζ = z \sqrt{ \\frac{Ω}{ν}}}$')
ax[1].set_xlabel('Dimensionless Velocity')


# In[3]:
'''
    Single Specie Transport to the rotating disk (no migration effects)
'''
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

''' Import Newman Figure 17.1 Data '''
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/LaTeX/Concentrated_Solution_Code/')
Newman_Data_Fig_17_1 = pd.read_table('Newman_Fig17_1_Data.txt', delimiter=',', header=0, skiprows=None)

ax[1].plot(Newman_Data_Fig_17_1['Position'], Newman_Data_Fig_17_1['Concentration'], 'ko', markerfacecolor='None')



'''
    Finite Volume
'''
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P5_DiluteSoln_RDE/FiniteVolume/SingleSpecie/Dimensionless')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data = data[data['Time'] == data['Time'].max()]

ax[1].plot(data['Position'], data['Conc'], linewidth=3, zorder=0, color = 'red', linestyle = '--')


'''
    Finite Difference
'''
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P5_DiluteSoln_RDE/FiniteDifference/SingleSpecie/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data = data[data['Time'] == data['Time'].max()]

ax[1].plot(data['Position'], data['Conc'], linewidth=6, zorder=-1, color = 'green')


os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/DiluteSoln/FiniteDifference/SingleSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['Conc'], 'b-', markerfacecolor='none')

ax[1].set_xlim(0,2.)
ax[1].set_ylim(0,1.05)

ax[1].set_xlabel('ξ')
ax[1].set_ylabel('Θ')





# In[4]:
'''
    Multi-Specie
        Migration effects included
'''

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P5_DiluteSoln_RDE/FiniteDifference/MultiSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['C_OH']*1e3, 'k-', linewidth=2)

ax_t.plot(data_['Position'], data_['C_Na']*1e3, 'k-', linewidth=2)
ax_t.plot(data_['Position'], data_['C_Cl']*1e3, 'k-', linewidth=2)
#

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P5_DiluteSoln_RDE/FiniteVolume/MultiSpecie/Dimensionless')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'ro', markerfacecolor='none')
ax[1].plot(data_['Position'], data_['C_OH']*1e3, 'ro', markerfacecolor='none')

ax_t.plot(data_['Position'], data_['C_Na']*1e3, 'ro', markerfacecolor='none')
ax_t.plot(data_['Position'], data_['C_Cl']*1e3, 'ro', markerfacecolor='none')




ax[1].set_xlabel('$\mathregular{ξ \ = \ y \ \left(\\frac{aν}{3D_R}\\right)^{1/3} \sqrt{ \\frac{Ω}{ν}} }$', fontsize = 20)
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0, 4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


ax[1].text(2, 0.000215, '$\mathregular{O_2}$', ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2, 0.000005, '$\mathregular{OH^-}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 1.00008, '$\mathregular{Na^+}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 0.9998, '$\mathregular{Cl^-}$', ha = 'left', va='bottom', fontsize = 15)



# In[4]:
'''
    Multi-Specie
        Migration effects included
'''

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/DiluteSoln/FiniteDifference/MultiSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['C_OH']*1e3, 'k-', linewidth=2)

ax_t.plot(data_['Position'], data_['C_Na']*1e3, 'k-', linewidth=2)
ax_t.plot(data_['Position'], data_['C_Cl']*1e3, 'k-', linewidth=2)
#




ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0, 4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


ax[1].text(2, 0.000215, '$\mathregular{O_2}$', ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2, 0.000005, '$\mathregular{OH^-}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 1.00008, '$\mathregular{Na^+}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 0.9998, '$\mathregular{Cl^-}$', ha = 'left', va='bottom', fontsize = 15)

# fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Dilute_Solution_Stagnant_Film_O2_Reduction.pdf', format='pdf', dpi=300, bbox_inches = "tight")


# In[4]:
'''
    Multi-Specie
        Migration effects included
'''

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/DiluteSoln/FiniteDifference/MultiSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['C_OH']*1e3, 'k-', linewidth=2)

ax_t.plot(data_['Position'], data_['C_Na']*1e3, 'k-', linewidth=2)
ax_t.plot(data_['Position'], data_['C_Cl']*1e3, 'k-', linewidth=2)


ax[2].plot(data_['Position'], data_['Phi']*1e6, 'k-', linewidth=2)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[2].set_ylabel('Φ (μV)')



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0, 4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


ax[1].text(2, 0.000215, '$\mathregular{O_2}$', ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2, 0.000005, '$\mathregular{OH^-}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 1.00008, '$\mathregular{Na^+}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 0.9998, '$\mathregular{Cl^-}$', ha = 'left', va='bottom', fontsize = 15)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Dilute_Solution_Stagnant_Film_O2_Reduction.pdf', format='pdf', dpi=300, bbox_inches = "tight")


# In[5]:
'''
    Stagnant Film Problem
        1 specie
        no migration
'''

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/DiluteSoln/FiniteDifference/SingleSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['Conc'], 'k-', markerfacecolor='none')

ax[1].set_xlim(0,2.)
ax[1].set_ylim(0,1.05)

ax[1].set_xlabel('ξ')
ax[1].set_ylabel('Θ')

# In[5]:
'''
    Stagnant Film Problem
        2 specie
        no migration
        Concentrated Solution Theory
'''

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
# ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_O2'], 'k-', markerfacecolor='none')
# ax_t.plot(data_['Position'], data_['x_H2O'], 'r-', markerfacecolor='none')


ax[2].plot(data_['Position'], data_['Flux_O2'], 'k-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_H2O'], 'k-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_O2'] + data_['Flux_H2O'], 'k--', markerfacecolor='none', linewidth=2)

ax[1].set_xlim(0,4.)
ax[1].set_ylim(ymin=0.0)
loc = matplotlib.ticker.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
ax[1].xaxis.set_major_locator(loc)

ax[1].set_xlabel('$\mathregular{η = \\frac{y}{ 2 \sqrt{\mathscr{D}_{0R} t} } }$', fontsize=20)
ax[1].set_ylabel('$\mathregular{x_{O2}}$')

ax[2].set_xlim(0,4.)
ax[2].xaxis.set_major_locator(loc)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{ 2 \sqrt{\mathscr{D}_{0R} t} } }$', fontsize=20)
ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)

ax[2].text(2.5, -1e-7,    '$\mathscr{N}_{\mathregular{O_2}}$', ha='left', va='top', fontsize=20)
ax[2].text(2.5, -4e-6, '$\mathscr{N}_{\mathregular{H_2O}}$', ha='left', va='bottom', fontsize=20)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Concentrated_Solution_Stagnant_Film_O2_in_H2O.pdf', format='pdf', dpi=300, bbox_inches = "tight")


# In[5]:
'''
    Stagnant Film Problem
        2 specie
        no migration
        Concentrated Solution Theory
'''

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
# ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_O2']*55.5, 'k-', markerfacecolor='none')
ax[1].plot(data_['Position'], data_['x_OH']*55.5, 'k-', markerfacecolor='none')

# ax_t = ax[1].twinx()
# ax_t.plot(data_['Position'], data_['x_Na'], 'k-', markerfacecolor='none')
# ax_t.plot(data_['Position'], data_['x_Cl'] + data_['Flux_H2O'], 'k--', markerfacecolor='none', linewidth=2)

ax[1].set_xlim(0,4.)
ax[1].set_ylim(ymin=0.0)

ax[1].set_xlabel('$\mathregular{ξ = \\frac{y}{ 2 \sqrt{\mathscr{D}_{0R} t} } }$', fontsize=20)
ax[1].set_ylabel('Mole Fraction')

ax[1].text(2, 0.000205, '$\mathregular{x_{O_2}}$', ha = 'left', va='bottom', fontsize = 20)
ax[1].text(2, 0.000045, '$\mathregular{x_{OH^-}}$', ha = 'left', va='bottom', fontsize = 20)

# ax[2].set_xlim(0,4.)
# ax[2].set_xlabel('ξ')
# ax[2].set_ylabel('$\mathscr{N}_i$')


# In[4]:
'''
    Multi-Specie
        Migration effects included
'''

rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_O2']*55.5, 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_OH']*55.5, 'k-', linewidth=2)

ax_t.plot(data_['Position'], data_['x_Na']*55.5, 'k-', linewidth=2)
ax_t.plot(data_['Position'], data_['x_Cl']*55.5, 'k-', linewidth=2)
#


ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Na'], 'k-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_H2O'], 'g-', markerfacecolor='none')
total_Flux = data_['Flux_O2'] + data_['Flux_OH'] + data_['Flux_Na'] + data_['Flux_Cl'] + data_['Flux_H2O']
# ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('$ \mathregular{[Na] \ and \ [Cl] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0, 4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


ax[1].text(2.9, 0.000211, '$\mathregular{O_2}$', ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2.9, 0.000005, '$\mathregular{OH}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 1.0+1e-6, '$\mathregular{Na}$', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(0.7, 1-2e-5, '$\mathregular{Cl}$', ha = 'left', va='top', fontsize = 15)


ax[2].text(1, -0.2e-5, '$\mathscr{N}_{\mathregular{O_2}}$', ha = 'left', va='top', fontsize = 20, color='red')
ax[2].text(1, 1.2e-5, '$\mathscr{N}_{\mathregular{OH}}$', ha = 'left', va='bottom', fontsize = 20, color='blue')
ax[2].text(3, 1.2e-5, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='bottom', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Concentrated_Solution_Stagnant_Film_O2_OH_NaCl_in_H2O_noMigration.pdf', format='pdf', dpi=300, bbox_inches = "tight")


# In[11]:
'''
    Concentrated Solution Theory
        Binary Salt
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/BinarySalt/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_Na'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_Cl'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_H2O'] - (55 - 2*1), 'k-', linewidth=2)
#
ax_t.plot(data_['Position'], data_['Phi']*1e3, '-', linewidth=2, color = 'red')

# ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Na'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_H2O'], 'g-', markerfacecolor='none')
total_Flux = data_['Flux_Na'] + data_['Flux_Cl'] + data_['Flux_H2O']
# ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('Concentration (mol/L)')

# ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('Solution Potential, $\mathregular{\Phi}$ (mV)', color='red')
ax_t.tick_params(axis='y', colors='red')
ax_t.spines['right'].set_color('red')
# ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0, 4.)
ax[1].set_ylim(0.45, 1.55)
ax_t.set_ylim(ymin=-100)


ax[1].text(2.5, 1, '$\mathregular{Na^+, \ Cl^-}$',  ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2.5, 0.51, '$\mathregular{H_2O}$ (+53)', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(3, -12, '$\mathregular{\Phi}$', ha = 'left', va='top', fontsize = 15, color='red')


ax[2].text(2.5, 1e-2, '$\mathscr{N}_{\mathregular{Cl^-}}$', ha = 'left', va='top', fontsize = 20, color='blue')
ax[2].text(2.5, -0.0075, '$\mathscr{N}_{\mathregular{Na^+}}$', ha = 'left', va='bottom', fontsize = 20, color='red')
ax[2].text(2.5, -0.022, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='bottom', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Concentrated_Solution_Stagnant_Film_BinarySalt_Migration.pdf', format='pdf', dpi=300, bbox_inches = "tight")

# In[11]:
'''
    Concentrated Solution Theory
        Binary Salt
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/BinarySalt/HighConc/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_Na'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_Cl'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_H2O'] - (55.5 - 2*10-5), 'k-', linewidth=2)
#
ax_t.plot(data_['Position'], data_['Phi']*1e3, '-', linewidth=2, color = 'red')

# ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Na'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_H2O'], 'g-', markerfacecolor='none')
total_Flux = data_['Flux_Na'] + data_['Flux_Cl'] + data_['Flux_H2O']
# ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('Concentration (mol/L)')

# ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('Solution Potential, $\mathregular{\Phi}$ (mV)', color='red')
ax_t.tick_params(axis='y', colors='red')
ax_t.spines['right'].set_color('red')
# ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0, 4.)
ax[1].set_ylim(4.5, 16)
ax_t.set_ylim(ymin=-100)


ax[1].text(2.5, 10.1, '$\mathregular{Na^+, \ Cl^-}$',  ha = 'left', va='bottom', fontsize = 15)
ax[1].text(2.5, 5.1, '$\mathregular{H_2O}$ (+25.5)', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(3, -12, '$\mathregular{\Phi}$', ha = 'left', va='top', fontsize = 15, color='red')


ax[2].text(2.5, 0.09, '$\mathscr{N}_{\mathregular{Cl^-}}$', ha = 'left', va='top', fontsize = 20, color='blue')
ax[2].text(2.5, -0.155, '$\mathscr{N}_{\mathregular{Na^+}}$', ha = 'left', va='bottom', fontsize = 20, color='red')
ax[2].text(2.5, -0.21, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='top', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Concentrated_Solution_Stagnant_Film_BinarySalt_HighConc_Migration.pdf', format='pdf', dpi=300, bbox_inches = "tight")


# In[11]:
'''
    Concentrated Solution Theory
        Binary Salt
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/BinarySalt/HighConc/SoluteInteractions/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_Na'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_Cl'], 'k-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_H2O'] - (55.5 - 2*10-5), 'k-', linewidth=2)
#
ax_t.plot(data_['Position'], data_['Phi']*1e3, '-', linewidth=2, color = 'red')

# ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Na'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_H2O'], 'g-', markerfacecolor='none')
total_Flux = data_['Flux_Na'] + data_['Flux_Cl'] + data_['Flux_H2O']
# ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('Concentration (mol/L)')

# ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('Solution Potential, $\mathregular{\Phi}$ (mV)', color='red')
ax_t.tick_params(axis='y', colors='red')
ax_t.spines['right'].set_color('red')
# ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0, 4.)
ax[1].set_ylim(4.5, 16)
ax_t.set_ylim(ymin=-400)


ax[1].text(5, 10.1, '$\mathregular{Na^+, \ Cl^-}$',  ha = 'left', va='bottom', fontsize = 15)
ax[1].text(5, 5.1, '$\mathregular{H_2O}$ (+25.5)', ha = 'left', va='bottom', fontsize = 15)
ax_t.text(8, -12, '$\mathregular{\Phi}$', ha = 'left', va='top', fontsize = 15, color='red')


ax[2].text(6, 0.2, '$\mathscr{N}_{\mathregular{Cl^-}}$', ha = 'left', va='top', fontsize = 20, color='blue')
ax[2].text(6, -0.33, '$\mathscr{N}_{\mathregular{Na^+}}$', ha = 'left', va='bottom', fontsize = 20, color='red')
ax[2].text(6, -0.44, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='top', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/')
fig.savefig('Concentrated_Solution_Stagnant_Film_BinarySalt_HighConc_SoluteInteractions_Migration.pdf', format='pdf', dpi=300, bbox_inches = "tight")



# In[11]:
'''
    Concentrated Solution Theory
        Two Binary Salts with Common Ion
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/CommonIon/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_Li'], 'r-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_BMP'], 'b-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_TFSI'], 'g-', linewidth=2)
#
ax_t.plot(data_['Position'], data_['Phi']*1e3, '-', linewidth=2, color = 'k')

# ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Li'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_BMP'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_TFSI'], 'g-', markerfacecolor='none')
total_Flux = data_['Flux_Li'] + data_['Flux_BMP'] + data_['Flux_TFSI']
ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('Concentration (mol/L)')

# ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('Solution Potential, $\mathregular{\Phi}$ (mV)', color='red')
ax_t.tick_params(axis='y', colors='red')
ax_t.spines['right'].set_color('red')
# ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0, 4.)
# ax[1].set_ylim(0.45, 1.55)
# ax_t.set_ylim(ymin=-0.4)


# ax[1].text(2.5, 1, '$\mathregular{Na^+, \ Cl^-}$',  ha = 'left', va='bottom', fontsize = 15)
# ax[1].text(2.5, 0.51, '$\mathregular{H_2O}$ (+53)', ha = 'left', va='bottom', fontsize = 15)
# ax_t.text(3, -12, '$\mathregular{\Phi}$', ha = 'left', va='top', fontsize = 15, color='red')
#
#
# ax[2].text(2.5, 1e-2, '$\mathscr{N}_{\mathregular{Cl^-}}$', ha = 'left', va='top', fontsize = 20, color='blue')
# ax[2].text(2.5, -0.0075, '$\mathscr{N}_{\mathregular{Na^+}}$', ha = 'left', va='bottom', fontsize = 20, color='red')
# ax[2].text(2.5, -0.022, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='bottom', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()

# In[11]:
'''
    Concentrated Solution Theory
        Two Binary Salts with Common Ion
'''
rows, columns = [1, 2]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/17_AppendixC/P6_InfiniteStagnantFilm/ConcSolnTheory/Dimensionless/Migration/CommonIon/x_0p5/')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_ = data[data['Time'] == data['Time'].max()]

ax[1].plot(data_['Position'], data_['x_Li'], 'r-', linewidth=2)
ax[1].plot(data_['Position'], data_['x_BMP'], 'b-', linewidth=2)
# ax[1].plot(data_['Position'], data_['x_TFSI']-0.4925, 'g-', linewidth=2)
#
ax_t.plot(data_['Position'], data_['Phi']*1e3, '-', linewidth=2, color = 'k')

# ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_Li'], 'r-', markerfacecolor='none')
ax[2].plot(data_['Position'], data_['Flux_BMP'], 'b-', markerfacecolor='none')
# ax[2].plot(data_['Position'], data_['Flux_TFSI'], 'g-', markerfacecolor='none')
# total_Flux = data_['Flux_Li'] + data_['Flux_BMP'] + data_['Flux_TFSI']
# ax[2].plot(data_['Position'], total_Flux, 'k--', markerfacecolor='none', linewidth=2)



ax[1].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)
ax[1].set_ylabel('Concentration (mol/L)', color='blue')
ax[1].tick_params(axis='y', colors='blue')
ax[1].spines['left'].set_color('blue')
# ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

ax_t.set_ylabel('Solution Potential, $\mathregular{\Phi}$ (mV)', color='black')
# ax_t.tick_params(axis='y', colors='red')
ax_t.spines['left'].set_color('blue')
# ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0, 4.)
# ax[1].set_ylim(0.45, 1.55)
# ax_t.set_ylim(ymin=-0.4)


# ax[1].text(2.5, 1, '$\mathregular{Na^+, \ Cl^-}$',  ha = 'left', va='bottom', fontsize = 15)
# ax[1].text(2.5, 0.51, '$\mathregular{H_2O}$ (+53)', ha = 'left', va='bottom', fontsize = 15)
# ax_t.text(3, -12, '$\mathregular{\Phi}$', ha = 'left', va='top', fontsize = 15, color='red')
#
#
# ax[2].text(2.5, 1e-2, '$\mathscr{N}_{\mathregular{Cl^-}}$', ha = 'left', va='top', fontsize = 20, color='blue')
# ax[2].text(2.5, -0.0075, '$\mathscr{N}_{\mathregular{Na^+}}$', ha = 'left', va='bottom', fontsize = 20, color='red')
# ax[2].text(2.5, -0.022, '$\mathscr{N}_{\mathregular{H_2O}}$', ha = 'left', va='bottom', fontsize = 20, color='green')

ax[2].set_ylabel('$\mathscr{N}_i$', fontsize=20)
ax[2].set_xlabel('$\mathregular{η = \\frac{y}{2 \sqrt{\mathscr{D}_{0R} t}} }$', fontsize=20)

fig.tight_layout()


# In[10]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)
ax_t = ax[1].twinx()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/FortranStuff/UpdateFortranMethods/Appendix_C/SteadyFillmat')

data = pd.read_table('Conc_profile_RDE_MultiComponent.txt', delim_whitespace=True, header=0, skiprows=[1])

ax[1].plot(data['Position'], data['conc_O2']*1e6)
ax[1].plot(data['Position'], data['conc_OH']*1e6)

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)


ax_t.plot(data['Position'], data['conc_Na']*1e3)
ax_t.plot(data['Position'], data['conc_Cl']*1e3)

ax_t.set_ylim(0.9992, 1.0002)


# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_AppendixC_P5_MultiSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data = data[data['Time'] == data['Time'].max()]
print(data.keys())

ax[1].plot(data['Xi'], data['C_O2']*1e6)
ax[1].plot(data['Xi'], data['C_OH']*1e6)
#
ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)
#
ax_t = ax[1].twinx()
ax_t.plot(data['Xi'], data['C_Na']*1e3)
ax_t.plot(data['Xi'], data['C_Cl']*1e3)
#
ax_t.set_ylim(0.9992, 1.0002)


# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/21_AppendixC_P5_FinDiff_MultSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data = data[data['Time'] == data['Time'].max()]
print(data.keys())

ax[1].plot(data['Position'], data['C_O2']*1e6)
ax[1].plot(data['Position'], data['C_OH']*1e6)
#
ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)
#
ax_t = ax[1].twinx()
ax_t.plot(data['Position'], data['C_Na']*1e3)
ax_t.plot(data['Position'], data['C_Cl']*1e3)
#
ax_t.set_ylim(0.9992, 1.0002)

# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/22_AppendixC_P6')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

ax[1].plot(data['Position'], data['C_O2']*1e3)
ax[1].plot(data['Position'], data['C_OH']*1e3)
#
ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)
#
ax_t = ax[1].twinx()
ax_t.plot(data['Position'], data['C_Na'])
ax_t.plot(data['Position'], data['C_Cl'])
#
ax_t.set_ylim(0.9992, 1.0002)

# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/22_AppendixC_P6')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])

ax[1].plot(data['Position'], data['C_H2O'])
# ax[1].plot(data['Position'], data['C_OH']*1e3)
#
ax[1].set_xlim(0,4.)
# ax[1].set_ylim(0,1)

# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/23_AppendixC_P6_MultiCompDiff')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax[1].plot(data['Position'], data['C_O2']*1e3)
ax[1].plot(data['Position'], data['C_OH']*1e3)
#
ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)
#
ax_t = ax[1].twinx()
ax_t.plot(data['Position'], data['C_Na'])
ax_t.plot(data['Position'], data['C_Cl'])
#
ax_t.set_ylim(0.9992, 1.0002)




# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/24_AppendixC_P6_DiluteSoln')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in [129]:
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'k-', linewidth=2)
    ax[1].plot(data_['Position'], data_['C_OH'], 'k-', linewidth=2)

    ax_t.plot(data_['Position'], data_['C_Na'], 'k-', linewidth=2)
    ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)
#

ax[1].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


# In[11]:
'''
    Compare profiles obtained through dilute solution theory to those obtained using concentrated solution theory
'''
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/24_AppendixC_P6_DiluteSoln')

data_dilute = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/22_AppendixC_P6')
data_conc = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
data_conc = data_conc[data_conc['Iter'] == data_conc['Iter'].max()]

ax_t = ax[1].twinx()

for it in [129]:
    data_ = data_dilute[data_dilute['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'ro', linewidth=2, markerfacecolor='none')
    ax[1].plot(data_['Position'], data_['C_OH'], 'ro', linewidth=2, markerfacecolor='none')

    ax_t.plot(data_['Position'], data_['C_Na'], 'ro', linewidth=2, markerfacecolor='none')
    ax_t.plot(data_['Position'], data_['C_Cl'], 'ro', linewidth=2, markerfacecolor='none')
#


ax[1].plot(data_conc['Position'], data_conc['C_O2'], 'k-', linewidth=4, zorder=-1)
ax[1].plot(data_conc['Position'], data_conc['C_OH'], 'k-', linewidth=4, zorder=-1)

ax_t.plot(data_conc['Position'], data_conc['C_Na'], 'k-', linewidth=4, zorder=-1)
ax_t.plot(data_conc['Position'], data_conc['C_Cl'], 'k-', linewidth=4, zorder=-1)




ax[1].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')
ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


# In[11]:
rows, columns = [3, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/22_AppendixC_P6')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in [129]:
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
    ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)

    ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
    ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)

ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0
# for it in [11]:
it = data['Iter'].max()
data_ = data[data['Iter'] == it]
ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)

ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)

i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)

ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)


ax[1].set_xticklabels([])
ax[2].set_xlim(0,4.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)

ax[2].text(1, 1.12e-5, '$\mathregular{OH^-}$', color='b')
ax[2].text(1, 0.4e-5, '$\mathregular{Cl^-}$', color='k')
ax[2].text(1.5, 0.05e-5, '$\mathregular{O_2}$', color='r')
ax[2].text(1.5, -0.7e-5, '$\mathregular{Na^+}$', color='g')
ax[2].text(3.8, 1.57e-5, '$\mathregular{-i / F = -\sum \ z_i \mathscr{N}_i }$', color='k', va='top', ha = 'right', fontsize=20)


ax[1].text(0.5, 8.1e-4, '$\mathregular{Na^+}$', color='g')
ax[1].text(0.5, 6.5e-4, '$\mathregular{Cl^-}$', color='k')

ax[1].text(0.5, 3.3e-4, '$\mathregular{OH^-}$', color='b')
ax[1].text(0.5, 0.5e-4, '$\mathregular{O_2}$', color='r')


# In[11]:
rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/24_AppendixC_P6_DiluteSoln')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in [0, 1, 3, 10, 129]:
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'r--', linewidth=1)
    ax[1].plot(data_['Position'], data_['C_OH']*1e3, 'b--', linewidth=1)

    ax_t.plot(data_['Position'], data_['C_Na'], 'g--', linewidth=1)
    ax_t.plot(data_['Position'], data_['C_Cl'], '--', color = 'purple', linewidth=1)
#

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1)
ax_t.set_ylim(0.9992, 1.0002)

# In[11]:
rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/24_AppendixC_P6_DiluteSoln')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in range(0,21):
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'r--', linewidth=1)
    ax[2].plot(data_['Position'], data_['C_OH']*1e3, 'b--', linewidth=1)

    ax[3].plot(data_['Position'], data_['C_Na'], 'g--', linewidth=1)
    ax[4].plot(data_['Position'], data_['C_Cl'], '--', color = 'purple', linewidth=1)
#

for i in range(1,5):
    ax[i].set_xlim(0,4.)

ax[1].set_ylim(0,1)
ax[2].set_ylim(0,1)

ax[3].set_ylim(0.9992, 1.0002)
ax[4].set_ylim(0.9992, 1.0002)

# In[11]:
rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/22_AppendixC_P6')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in range(0,21):
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2']*1e3, 'r--', linewidth=1)
    ax[2].plot(data_['Position'], data_['C_OH']*1e3, 'b--', linewidth=1)

    ax[3].plot(data_['Position'], data_['C_Na'], 'g--', linewidth=1)
    ax[4].plot(data_['Position'], data_['C_Cl'], '--', color = 'purple', linewidth=1)
#

for i in range(1,5):
    ax[i].set_xlim(0,4.)

ax[1].set_ylim(0,1)
ax[2].set_ylim(0,1)

ax[3].set_ylim(0.9992, 1.0002)
ax[4].set_ylim(0.9992, 1.0002)

# In[11]:
rows, columns = [3, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/25_AppendixC_P6_DiluteSoln_with_Fluxes')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

for it in [129]:
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
    ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)

    ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
    ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)

ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)


z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0
# for it in [11]:
it = data['Iter'].max()
data_ = data[data['Iter'] == it]
ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)

ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)

i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)

ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)


ax[1].set_xticklabels([])
ax[2].set_xlim(0,4.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)

ax[2].text(1, 1.12e-5, '$\mathregular{OH^-}$', color='b')
ax[2].text(1, 0.4e-5, '$\mathregular{Cl^-}$', color='k')
ax[2].text(1.5, 0.05e-5, '$\mathregular{O_2}$', color='r')
ax[2].text(1.5, -0.7e-5, '$\mathregular{Na^+}$', color='g')
ax[2].text(3.8, 1.57e-5, '$\mathregular{-i / F = -\sum \ z_i \mathscr{N}_i }$', color='k', va='top', ha = 'right', fontsize=20)


ax[1].text(0.5, 8.1e-4, '$\mathregular{Na^+}$', color='g')
ax[1].text(0.5, 6.5e-4, '$\mathregular{Cl^-}$', color='k')

ax[1].text(0.5, 3.3e-4, '$\mathregular{OH^-}$', color='b')
ax[1].text(0.5, 0.5e-4, '$\mathregular{O_2}$', color='r')

# fig.tight_layout()


# In[11]:
rows, columns = [3, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/25_AppendixC_P6_DiluteSoln_with_Fluxes')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0


for it in range(0,2):
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
    ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)

    ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
    ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)




    ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
    ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)

    ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
    ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)


    # i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

    # ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)
    #
    ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)





ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,4.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)



ax[1].set_xticklabels([])
ax[2].set_xlim(0,4.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)

# ax[2].text(1, 1.12e-5, '$\mathregular{OH^-}$', color='b')
# ax[2].text(1, 0.4e-5, '$\mathregular{Cl^-}$', color='k')
# ax[2].text(1.5, 0.05e-5, '$\mathregular{O_2}$', color='r')
# ax[2].text(1.5, -0.7e-5, '$\mathregular{Na^+}$', color='g')
# ax[2].text(3.8, 1.57e-5, '$\mathregular{-i / F = -\sum \ z_i \mathscr{N}_i }$', color='k', va='top', ha = 'right', fontsize=20)
#
#
# ax[1].text(0.5, 8.1e-4, '$\mathregular{Na^+}$', color='g')
# ax[1].text(0.5, 6.5e-4, '$\mathregular{Cl^-}$', color='k')
#
# ax[1].text(0.5, 3.3e-4, '$\mathregular{OH^-}$', color='b')
# ax[1].text(0.5, 0.5e-4, '$\mathregular{O_2}$', color='r')

# fig.tight_layout()



# In[11]:
rows, columns = [3, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_AppendixC_P6_DiluteSoln_with_Fluxes2')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0


for it in range(0,2):
    data_ = data[data['Iter'] == it]
    ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
    ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)

    ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
    ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)




    ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
    ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)

    ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
    ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)


    i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

    ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)
    #
    ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)





ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,6.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)



ax[1].set_xticklabels([])
ax[2].set_xlim(0,6.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)

# ax[2].text(1, 1.12e-5, '$\mathregular{OH^-}$', color='b')
# ax[2].text(1, 0.4e-5, '$\mathregular{Cl^-}$', color='k')
# ax[2].text(1.5, 0.05e-5, '$\mathregular{O_2}$', color='r')
# ax[2].text(1.5, -0.7e-5, '$\mathregular{Na^+}$', color='g')
# ax[2].text(3.8, 1.57e-5, '$\mathregular{-i / F = -\sum \ z_i \mathscr{N}_i }$', color='k', va='top', ha = 'right', fontsize=20)
#
#
# ax[1].text(0.5, 8.1e-4, '$\mathregular{Na^+}$', color='g')
# ax[1].text(0.5, 6.5e-4, '$\mathregular{Cl^-}$', color='k')
#
# ax[1].text(0.5, 3.3e-4, '$\mathregular{OH^-}$', color='b')
# ax[1].text(0.5, 0.5e-4, '$\mathregular{O_2}$', color='r')

# fig.tight_layout()


# In[11]:
rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/27_AppendixC_P6_DiluteSoln_OneSpecie')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0


it = data['Iter'].max()
data_ = data[data['Iter'] == it]
ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
# ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)
#
# ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
# ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)




ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)
#
# ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
# ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)


# i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

# ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)
    #
    # ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)





ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

# ax[1].set_xlim(0,6.)
ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)



ax[1].set_xticklabels([])
# ax[2].set_xlim(0,6.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)


# In[11]:
rows, columns = [2, 1]
ax, fig = axes(1, rows, columns)

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/30_AppendixC_P6_DiluteSoln_OneSpecie_Dual3')

data = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=0, skiprows=[1])
# data = data[data['Time'] == data['Time'].max()]

ax_t = ax[1].twinx()

z_O2 = 0.0
z_OH = -1.0
z_Cl = -1.0
z_Na = +1.0


it = data['Iter'].max()
it = 2
data_ = data[data['Iter'] == it]
ax[1].plot(data_['Position'], data_['C_O2'], 'r-', linewidth=2)
# ax[1].plot(data_['Position'], data_['C_OH'], 'b-', linewidth=2)
#
# ax_t.plot(data_['Position'], data_['C_Na'], 'g-', linewidth=2)
# ax_t.plot(data_['Position'], data_['C_Cl'], 'k-', linewidth=2)




ax[2].plot(data_['Position'], data_['Flux_O2'], 'r-', linewidth=2)
# ax[2].plot(data_['Position'], data_['Flux_OH'], 'b-', linewidth=2)
#
# ax[2].plot(data_['Position'], data_['Flux_Na'], 'g-', linewidth=2)
# ax[2].plot(data_['Position'], data_['Flux_Cl'], 'k-', linewidth=2)


# i_F = data_['Flux_O2']*z_O2 + data_['Flux_OH']*z_OH + data_['Flux_Na']*z_Na + data_['Flux_Cl']*z_Cl

# ax[2].plot(data_['Position'], -i_F, 'k--', linewidth=2)
    #
    # ax[3].plot(data_['Position'], data_['Phi'], 'k--', linewidth=2)





ax[1].set_ylabel('$ \mathregular{[O_2] \ and \ [OH^-] }$')

ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# ax[1].yaxis.major.formatter.set_powerlimits((-3,-3))

ax_t.set_ylabel('$ \mathregular{[Na^+] \ and \ [Cl^-] }$')
ax_t.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax[1].set_xlim(0,6.)
# ax[1].set_ylim(0,1e-3)
ax_t.set_ylim(0.9992, 1.0002)



ax[1].set_xticklabels([])
ax[2].set_xlim(0,6.)
ax[2].set_xlabel('$\mathregular{η = y / (2\sqrt{D_{R}t} ) }$')

ax[2].set_ylabel('$ \mathscr{N}_i $', fontsize=20)





# In[1000]:
