import sys
sys.path.insert(0, '/Users/nicholasbrady/Documents/School/Academic/West Research/Projects/')
from PythonHelperFunctions import *
import RedlichKister as RKA
# from RedlichKister import *
plot_parameters()
csfont = {'fontname':'Serif'}
# latfont = {'Computer Modern Roman'}
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# matplotlib.rc('text', usetex=False)
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# matplotlib.rc('text', usetex=True)
# https://matplotlib.org/users/customizing.html
# matplotlib.use('cairo')
matplotlib.rc('xtick', top=True, bottom=True, direction='in')
matplotlib.rc('ytick', left=True, right=True, direction='in')
plt.rc("axes.spines", top=True, right=True)
matplotlib.rc('axes', edgecolor='k')
# matplotlib.rc('text', usetex=False)
# matplotlib.rc('font',**{'family':'serif', 'serif':['Times']})
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

import impedance
from impedance import preprocessing
from impedance.visualization import plot_nyquist
from impedance.models.circuits import CustomCircuit
from tqdm.notebook import tqdm as tqdm


def plot_Volt_Current(data, ax, title):

    ax.set_title(title)

    ax.plot(data['time/s']/3600, data['Ewe/V']*1e3, color='b', linewidth=2)
    ax.tick_params(axis='y', colors='blue')

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Voltage (mV)', color='blue')

    ax_t = ax.twinx()
    ax_t.plot(data['time/s']/3600, data['control/mA'] * 1e3, 'r', linewidth = 1)
    ax_t.set_ylabel('Current (μΑ)', color = 'r')
    ax_t.tick_params(axis='y', colors='red')
    ax_t.spines['right'].set_color('red')
    ax_t.spines['left'].set_color('blue')

    fig.tight_layout()

    return ax_t


color_list = ['blue', 'brown', 'darkorange', 'green', 'purple', 'black']

Fconst = 96485
Rigc   = 8.314
Temp   = 298

# In[1]:
LiPF6_in_EC_EMC_dir = '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_LiPF6_in_EC_EMC/'
os.chdir(LiPF6_in_EC_EMC_dir)
os.chdir('No_Current/')

data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1], index_col=None)
print(data.keys())
time = np.unique(data['Time'])

ax, fig = axes(1, rows=2, columns=3)

end_of_discharge = np.unique(data['Time'][data['State'] != 'R'])[-1]
rest_ind = np.unique(data['Time'][data['State'] == 'R'])[0]

for t in [0.0, 0.1, 1.0, 5, 10, end_of_discharge]:
    ind = data['Time'] == t
    ax[1].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[2].plot(data['Position'][ind],
                data['velocity'][ind])

    ax[3].plot(data['Position'][ind],
                data['Density'][ind])


ax[1].set_ylabel('Concentration [M]')
ax[1].set_xlabel('Position [cm]')

ax[2].set_ylabel('Mass-Average Velocity [cm/s]')
ax[2].set_xlabel('Position [cm]')

ax[3].set_ylabel('Density [g/mL]')
ax[3].set_xlabel('Position [cm]')



for t in [end_of_discharge+1, end_of_discharge+5, end_of_discharge+10, end_of_discharge+20]:
    ind = data['Time'] == t
    ax[4].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[5].plot(data['Position'][ind],
                data['velocity'][ind])

    ax[6].plot(data['Position'][ind],
                data['Density'][ind])


ax[4].set_ylabel('Concentration [M]')
ax[4].set_xlabel('Position [cm]')

ax[5].set_ylabel('Mass-Average Velocity [cm/s]')
ax[5].set_xlabel('Position [cm]')

ax[6].set_ylim(ax[3].get_ylim())
ax[6].set_ylabel('Density [g/mL]')
ax[6].set_xlabel('Position [cm]')

# Add text to plot indicating direction of flux (Current ON) and time
ax[1].text(0.5, 2.0, '$𝐍$', ha = 'left', va='center', color='k', fontsize=20)
ax[1].annotate("", xy=(0.7, 2.0), xytext=(0.55, 2.0), arrowprops=dict(arrowstyle="->"))
ax[1].annotate("time", xy=(0.1, 1.01), xytext=(0.3, 1.75), arrowprops=dict(arrowstyle="<-"))

# Add text to plot indicating zero Flux (Current OFF) and direction of time
ax[4].text(0.5, 1.8, '$𝐍 = 0$', ha = 'left', va='center', color='k', fontsize=20)
ax[4].annotate("time", xy=(0.1, 1.0), xytext=(0.25, 1.7), arrowprops=dict(arrowstyle="->"))

fig.savefig('Conc_Profiles_ConcSoln_NoCurrent.pdf', format='pdf', dpi=100, bbox_inches = "tight")

# In[]:
'''
    In and effort to diagnose some problems in the code, I first made a concentrated solution theory code without electrochemical effects - just diffusion and convection (ρ𝐯)
    Equations:                                  Boundary West                   Boundary East
    (1) ∂cₑ/∂t = ∂(ρωₑ/MWₑ)/∂t                  N₁ = 𝐢/F                         N₁ = 𝐢/F
    (2) ∂ρ/∂t = -∇⋅(ρ𝐯)                         𝐍₂ = 𝐍₁ * MWₑ                   𝐍₂ = ρ𝐯

    The top row plots the two dependent variables (cₑ, 𝐯) when the "current" is ON
    The bottom rows plots the two dependent variables (cₑ, 𝐯) when the "current" is OFF
    The concentration profiles (left) exhibit normal behavior - polarization when the current is ON and relaxation when the current is OFF

    The mass-averaged velocity profiles are a bit less intuitive (at least for me). If we reach a steady-state, then the mass-average velocity will be a constant value (MWₑ⋅𝐢/F), i.e. just the flowrate of the salt. But in the transient we see that the mass-average velocity is actually less than the steady-state value. If we look at the evolution in time of the density profiles (right), we can see that the solution is becoming more dense at the left boundary and less dense at the right boundary. So the mass-average velocity is less than the mass flow of MWₑ⋅𝐢/F because mass needs to flow to the left to densify the solution on the left and de-density the solution on the right.
'''

# In[2]:
os.chdir(LiPF6_in_EC_EMC_dir)
os.chdir('With_Current/')

data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1], index_col=None)
print(data.keys())
time = np.unique(data['Time'])

ax, fig = axes(1, rows=2, columns=3)
ax2, fig2 = axes(2, rows=1, columns=1)

end_of_discharge = np.unique(data['Time'][data['State'] != 'R'])[-1]
rest_ind = np.unique(data['Time'][data['State'] == 'R'])[0]

for t in [0.0, 0.1, 1.0, 5, 10, end_of_discharge]:
    ind = data['Time'] == t
    ax[1].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[2].plot(data['Position'][ind][:],
                data['velocity'][ind][:])

    ax[3].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)

ax[1].set_ylabel('Concentration [M]')
ax[1].set_xlabel('Position [cm]')

ax[2].set_ylabel('Mass-Average Velocity [cm/s]')
ax[2].set_xlabel('Position [cm]')

ax[3].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[3].set_xlabel('Position [cm]')


for t in [end_of_discharge+1, end_of_discharge+5, end_of_discharge+10, end_of_discharge+20]:
    ind = data['Time'] == t
    ax[4].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[5].plot(data['Position'][ind][:],
                data['velocity'][ind][:])

    ax[6].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)

ax[4].set_ylabel('Concentration [M]')
ax[4].set_xlabel('Position [cm]')

ax[5].set_ylabel('Mass-Average Velocity [cm/s]')
ax[5].set_xlabel('Position [cm]')

ax[6].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[6].set_xlabel('Position [cm]')


ind = data['Time'] == rest_ind
rest_c_LiPF6 = data['c_LiPF6'][ind]*1e3
ax2[1].plot(data['c_LiPF6'][ind]*1e3, data['Φ_2'][ind]*1e3, 'ko', fillstyle='none')

print(data['c_LiPF6'][ind].values[6]*1e3, data['Φ_2'][ind].values[6]*1e3)

RKAs = [ -68.7126396712847,
         -434.1883076829836 ,
         -225.62937774695456,
         -227.21395278116688,
         -218.01369472872466,
         -150.53946793731302]

x_hat  = np.linspace(1e-3, 3, 50)
V_hat = RKA.Volt_RKA_OCP(x_hat, RKAs, cimax=5)
offset = RKA.Volt_RKA_OCP(min(rest_c_LiPF6), RKAs, cimax=5)

ax2[1].plot(x_hat, V_hat - offset, 'r--', zorder=-100, linewidth=1)
ax2[1].set_ylabel('$\mathregular{\Phi_2}$ (mV)')
ax2[1].set_xlabel('$\mathregular{c_{LiPF_6}}$ (mol/L)')

ax2[1].text(.4, 40, 'Fit', ha = 'left', va='center', color = 'red', fontsize=20)
ax2[1].plot([0.15, 0.38], [15, 35], 'r-', linewidth=1)

ax2[1].text(1.7, -100, 'Simulated Data', ha = 'left', va='bottom', color = 'black', fontsize=20)
ax2[1].plot([data['c_LiPF6'][ind].values[6]*1e3, 1.68], [data['Φ_2'][ind].values[6]*1e3, -90],
                'k-', linewidth=1)

ax[1].text(0.5, 3.0, '$\mathregular{i_{app}}$', ha = 'left', va='center', color='k', fontsize=20)
ax[1].annotate("", xy=(0.7, 3.0), xytext=(0.6, 3.0), arrowprops=dict(arrowstyle="->"))
ax[1].annotate("time", xy=(0.1, 1.01), xytext=(0.3, 1.75), arrowprops=dict(arrowstyle="<-"))

ax[4].text(0.5, 2.3, '$\mathregular{i_{app} = 0}$', ha = 'left', va='center', color='k', fontsize=20)
ax[4].annotate("time", xy=(0.1, 1.0), xytext=(0.25, 1.7), arrowprops=dict(arrowstyle="->"))

fig.tight_layout()
fig2.tight_layout()

fig.savefig('Conc_Profiles_ConcSoln_WithCurrent.pdf', format='pdf', dpi=100, bbox_inches = "tight")
fig2.savefig('Phi_vs_Conc_Conc_Soln.pdf', format='pdf', dpi=100, bbox_inches = "tight")
#

# In[3]
dir_DiluteSoln = LiPF6_in_EC_EMC_dir + 'Dilute_Solution_Theory/'
os.chdir(dir_DiluteSoln + 't_Li_LessThan_Half/')

data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1], index_col=None)
print(data.keys())
time = np.unique(data['Time'])

ax, fig   = axes(1, rows=2, columns=2)
ax2, fig2 = axes(2, rows=2, columns=2)
ax3, fig3 = axes(3, rows=1, columns=2)

end_of_discharge = np.unique(data['Time'][data['State'] != 'R'])[-1]
rest_ind = np.unique(data['Time'][data['State'] == 'R'])[0]


for t in [0.0, 0.1, 1.0, 5, 10, end_of_discharge]:
    ind = data['Time'] == t
    ax[1].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[2].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax[1].set_ylabel('Concentration [M]')
ax[1].set_xlabel('Position [cm]')


ax[2].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[2].set_xlabel('Position [cm]')


for t in [end_of_discharge+1, end_of_discharge+5, end_of_discharge+10, end_of_discharge+20]:
    ind = data['Time'] == t
    ax[3].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[4].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax[3].set_ylabel('Concentration [M]')
ax[3].set_xlabel('Position [cm]')


ax[4].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[4].set_xlabel('Position [cm]')


ax[1].text(0.5, 2.3, '$\mathregular{i_{app}}$', ha = 'left', va='top', color='k', fontsize=20)
ax[1].annotate("", xy=(0.7, 2.2), xytext=(0.6, 2.2), arrowprops=dict(arrowstyle="->"))
ax[1].text(0.5, 2.1, '$\mathregular{t_{Li^+} < 0.5}$', ha = 'left', va='top', color='k', fontsize=20)

ax[3].text(0.5, 1.95, '$\mathregular{i_{app} = 0}$', ha = 'left', va='top', color='k', fontsize=20)
ax[3].text(0.5, 1.8, '$\mathregular{t_{Li^+} < 0.5}$', ha = 'left', va='top', color='k', fontsize=20)


ind = data['Time'] == rest_ind
rest_c_LiPF6 = data['c_LiPF6'][ind]*1e3
ax3[1].plot(data['c_LiPF6'][ind]*1e3, data['Φ_2'][ind]*1e3, 'ko-', linewidth=1)
ax3[1].set_ylabel('$\mathregular{\Phi_2}$ (mV)')
ax3[1].set_xlabel('$\mathregular{c_{LiPF_6}}$ (mol/L)')


os.chdir(LiPF6_in_EC_EMC_dir)
os.chdir('Dilute_Solution_Theory/t_Li_GreaterThan_Half')

data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1], index_col=None)
print(data.keys())
time = np.unique(data['Time'])

end_of_discharge = np.unique(data['Time'][data['State'] != 'R'])[-1]
rest_ind = np.unique(data['Time'][data['State'] == 'R'])[0]


for t in [0.0, 0.1, 1.0, 5, 10, end_of_discharge]:
    ind = data['Time'] == t
    ax2[1].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax2[2].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax2[1].set_ylabel('Concentration [M]')
ax2[1].set_xlabel('Position [cm]')


ax2[2].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax2[2].set_xlabel('Position [cm]')


for t in [end_of_discharge+1, end_of_discharge+5, end_of_discharge+10, end_of_discharge+20]:
    ind = data['Time'] == t
    ax2[3].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax2[4].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax2[3].set_ylabel('Concentration [M]')
ax2[3].set_xlabel('Position [cm]')


ax2[4].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax2[4].set_xlabel('Position [cm]')

ax2[1].text(0.5, 1.45, '$\mathregular{i_{app}}$', ha = 'left', va='top', color='k', fontsize=20)
ax2[1].annotate("", xy=(0.7, 1.42), xytext=(0.6, 1.42), arrowprops=dict(arrowstyle="->"))
ax2[1].text(0.5, 1.35, '$\mathregular{t_{Li^+} > 0.5}$', ha = 'left', va='top', color='k', fontsize=20)



ax2[3].text(0.5, 1.35, '$\mathregular{i_{app} = 0}$', ha = 'left', va='top', color='k', fontsize=20)
ax2[3].text(0.5, 1.28, '$\mathregular{t_{Li^+} > 0.5}$', ha = 'left', va='top', color='k', fontsize=20)

ind = data['Time'] == rest_ind
rest_c_LiPF6 = data['c_LiPF6'][ind]*1e3
ax3[2].plot(data['c_LiPF6'][ind]*1e3, data['Φ_2'][ind]*1e3, 'rs--', linewidth=1)
ax3[2].set_ylabel('$\mathregular{\Phi_2}$ (mV)')
ax3[2].set_xlabel('$\mathregular{c_{LiPF_6}}$ (mol/L)')

ax3[1].text(0.4, 17, '$\mathregular{t_{Li^+} < 0.5}$', ha='left', va='top', color='k', fontsize = 20)
ax3[1].text(0.4, 15, '$\mathregular{D_{Li^+} < \ D_{PF_6^-}}$', ha='left', va='top', color='k', fontsize = 20)
ax3[1].text(0.4, 13, '$\mathregular{d\Phi/dc > 0}$', ha='left', va='top', color='k', fontsize = 20)

ax3[2].text(1.15, 0, '$\mathregular{t_{Li^+} > 0.5}$', ha='left', va='top', color='r', fontsize = 20)
ax3[2].text(1.15, -0.75, '$\mathregular{D_{Li^+} > \ D_{PF_6^-}}$', ha='left', va='top', color='r', fontsize = 20)
ax3[2].text(1.15, -0.75*2, '$\mathregular{d\Phi/dc < 0}$', ha='left', va='top', color='r', fontsize = 20)

fig.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

ax[1].annotate("time", xy=(0.1, 1.01), xytext=(0.3, 1.75), arrowprops=dict(arrowstyle="<-"))
ax[3].annotate("time", xy=(0.1, 1.0), xytext=(0.25, 1.7), arrowprops=dict(arrowstyle="->"))
ax2[1].annotate("time", xy=(0.1, 1.01), xytext=(0.3, 1.3), arrowprops=dict(arrowstyle="<-"))
ax2[3].annotate("time", xy=(0.1, 1.01), xytext=(0.25, 1.25), arrowprops=dict(arrowstyle="->"))

os.chdir(dir_DiluteSoln)
fig.savefig('Conc_Profiles_t_Li_LessThan_Half.pdf', format='pdf', dpi=100, bbox_inches = "tight")
fig2.savefig('Conc_Profiles_t_Li_GreaterThan_Half.pdf', format='pdf', dpi=100, bbox_inches = "tight")
fig3.savefig('Phi_vs_Conc_Profiles.pdf', format='pdf', dpi=100, bbox_inches = "tight")

#

# In[4]:
os.chdir(dir_DiluteSoln)

data = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=0, skiprows=[1], index_col=None)
print(data.keys())
time = np.unique(data['Time'])

ax, fig   = axes(1, rows=2, columns=2)

end_of_discharge = np.unique(data['Time'][data['State'] != 'R'])[-1]
rest_ind = np.unique(data['Time'][data['State'] == 'R'])[0]


for t in [0.0, 0.1, 1.0, 5, 10, end_of_discharge]:
    ind = data['Time'] == t
    ax[1].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[2].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax[1].set_ylabel('Concentration [M]')
ax[1].set_xlabel('Position [cm]')


ax[2].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[2].set_xlabel('Position [cm]')


for t in [end_of_discharge+1, end_of_discharge+5, end_of_discharge+10, end_of_discharge+20]:
    ind = data['Time'] == t
    ax[3].plot(data['Position'][ind],
                data['c_LiPF6'][ind]*1e3)

    ax[4].plot(data['Position'][ind],
                data['Φ_2'][ind]*1e3)
#
ax[3].set_ylabel('Concentration [M]')
ax[3].set_xlabel('Position [cm]')


ax[4].set_ylabel('$\mathregular{\Delta \Phi_2}$ (mV)')
ax[4].set_xlabel('Position [cm]')

ax[1].text(0.5, 2.5, '$\mathregular{i_{app}}$', ha = 'left', va='center', color='k', fontsize=20)
ax[1].annotate("", xy=(0.7, 2.5), xytext=(0.6, 2.5), arrowprops=dict(arrowstyle="->"))
ax[1].annotate("time", xy=(0.1, 1.01), xytext=(0.25, 2.0), arrowprops=dict(arrowstyle="<-"))

ax[3].text(0.5, 2.1, '$\mathregular{i_{app} = 0}$', ha = 'left', va='center', color='k', fontsize=20)
ax[3].annotate("time", xy=(0.1, 1.0), xytext=(0.25, 1.7), arrowprops=dict(arrowstyle="->"))

fig.tight_layout()
fig.savefig('Conc_Profiles_t_Li_Measured.pdf', format='pdf', dpi=100, bbox_inches = "tight")
