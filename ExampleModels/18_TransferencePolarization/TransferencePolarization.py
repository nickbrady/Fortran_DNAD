# In[0]:
# import libraries and helper functions
# Please see PythonHelperFunctions for what exactly is imported
from cmath import polar
import sys
sys.path.insert(0, '/Users/nicholasbrady/Documents/School/Academic/West Research/Projects/')
from PythonHelperFunctions import *
plot_parameters()
csfont = {'fontname':'Serif'}
# latfont = {'Computer Modern Roman'}
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

def read_Fortran(file):
    return pd.read_table(file, delim_whitespace=True, header=[0], skiprows=[1])


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

# In[1]:
# Analytical solution 
# 𝜕𝑐/𝜕𝑡=𝐷 (𝜕²𝑐)/(𝜕𝑥² )
# 𝑐(𝑡=0,𝑥) = 𝑐_∞
# 𝑐(𝑡,𝑥=∞) = 𝑐_∞
# 𝜕𝑐/𝜕𝑥(𝑥=0) = 𝐾[𝑢(𝑡)−𝑢(𝑡−𝑡_𝑖 )]
# 𝐾 = -(1 - t_Li)I/(z_Li nu_Li F D)

rows, columns = [2, 2]
ax, fig = axes(1, rows=rows, columns=columns)

def func(K, D, t, x, ti=np.inf):
    '''
        https://math.stackexchange.com/questions/4213141/inverse-laplace-transform-of-exponential-function/4217293#4217293
        Δc = c(t,x) - c_inf =          -K { \sqrt(4Dt/π)      * exp(-x^2/(4Dt))      - x erfc(x/sqrt(4Dt)) - 
                                u(t - ti)*[ \sqrt(4D(t-ti)/π) * exp(-x^2/(4D(t-ti))) - x erfc(x/sqrt(4D(t-ti)))] }

        u(t-ti) = 0 if t <= ti
                = 1 if t >  ti
    '''
    def f(D, t, x):
        return np.sqrt(4*D*t/np.pi) * np.exp(-x**2/(4*D*t)) - x * scipy.special.erfc(x/np.sqrt(4*D*t))

    f_less_ti = f(D, t, x)

    if t <= ti:
        return -K * f_less_ti

    else: # t > ti
                    # f_lt_ti   - f_gt_ti
        return -K * (f_less_ti - f(D, t - ti, x))


def func_x0(K, D, t, ti=np.inf):
    '''
        Δc = c(t,x=0) - c_inf = -K { \sqrt(4Dt/π) - u(t - ti)*[ \sqrt(4D(t-ti)/π) ] }
    '''
    def f(D, t):
        return np.sqrt(4*D*t/np.pi)

    delta_c = np.where(t <= ti, f(D, t), f(D, t) - f(D, t-ti))
    # delta_c = np.where(t <= ti, np.sqrt(4*D*t/np.pi), np.sqrt(4*D*t/np.pi) - np.sqrt(4*D*(t-ti)/np.pi))

    return -K * delta_c


K = -1
D = 1

def makeFigure():

    t_ = np.logspace(-1, 1, 4)

    c_inf = 1.0

    x = np.linspace(0, 10, 200)

    t_ = [1e-2, 0.25, 1, 4, 10]
    for t in t_:
        delta_C = func(K, D, t, x)
        ax[1].plot(x, delta_C)


    t_ = np.logspace(-2, 1)
    for t in t_:
        delta_C = func(K, D, t, x)
        ax[3].plot(t, delta_C[0], 'ko')

    ax[3].plot(t_, func_x0(K, D, t_), 'r')


    t_ = [10, 11, 15, 30, 100, 1e5]
    for t in t_:
        delta_C = func(K, D, t, x, ti=10)
        ax[2].plot(x, delta_C)

    t_ = np.logspace(1, 2)
    ti = 10
    for t in t_:
        delta_C = func(K, D, t, x, ti=ti)
        ax[4].plot(t, delta_C[0], 'ko')

    delta_C = func_x0(K, D, t_, ti)
    ax[4].plot(t_, delta_C, 'r')

    ylims = ax[3].get_ylim()
    ax[4].set_ylim(ylims)

    ax[1].set_xlim(xmin=-1.9)
    ax[2].set_xlim(xmin=-1.9)

    _=ax[2].set_yticklabels([])
    _=ax[4].set_yticklabels([])

    ax[1].set_ylabel('$\mathregular{\Delta c}$', fontsize=20)

    ax[1].set_xlabel('$\mathregular{x}$', fontsize=20)
    ax[2].set_xlabel('$\mathregular{x}$', fontsize=20)

    ax[1].set_title('Current ON', y=0.9, fontsize=20, fontweight='bold')
    ax[2].set_title('Current OFF', y=0.9, fontsize=20, fontweight='bold')
    ax[1].text(10, 3.5, '$t \leq t_i$', va='top', ha='right', fontsize=20)
    ax[2].text(10, 3.5, '$t > t_i$', va='top', ha='right', fontsize=20)


    ax[3].set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=20)

    ax[3].set_xlabel('$t$', fontsize=20)
    ax[4].set_xlabel('$t$', fontsize=20)

    ax[3].text(0,  3.5, '$t \leq t_i$', va='top', ha='left', fontsize=20)
    ax[4].text(15, 3.5, '$t > t_i$', va='top', ha='left', fontsize=20)

    ax[1].text(-0.1, 0,    '$t = 0.01$', va='bottom', ha='right', fontsize=13, color='blue')
    ax[1].text(-0.1, .5,   '$t = 0.25$', va='bottom', ha='right', fontsize=13, color='darkorange')
    ax[1].text(-0.1, 1.1,  '$t = 1$',    va='bottom', ha='right', fontsize=13, color='green')
    ax[1].text(-0.1, 2.2,  '$t = 4$',    va='bottom', ha='right', fontsize=13, color='red')
    ax[1].text(-0.1, 3.5,  '$t = 10$',   va='bottom', ha='right', fontsize=13, color='purple')

    ax[2].text(-0.1, 3.5, '$t = 10$',     va='bottom', ha='right', fontsize=13, color='blue')
    ax[2].text(-0.1, 2.6, '$t = 11$',     va='bottom', ha='right', fontsize=13, color='darkorange')
    ax[2].text(-0.1, 1.8, '$t = 15$',     va='bottom', ha='right', fontsize=13, color='green')
    ax[2].text(-0.1, 1.1, '$t = 30$',     va='bottom', ha='right', fontsize=13, color='red')
    ax[2].text(-0.1, 0.5, '$t = 100$',    va='bottom', ha='right', fontsize=13, color='purple')
    ax[2].text(-0.1, 0,   '$t = \infty$', va='center', ha='right', fontsize=13, color='brown')

    ax[1].set_yticks([0, 1, 2, 3])
    ax[1].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.5))
    ax[1].tick_params(which='minor', length=4)

    ax[3].set_yticks([0, 1, 2, 3])
    ax[3].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.5))
    ax[3].tick_params(which='minor', length=4)

makeFigure()

ax[3].text(0, 2.5, '$Δc = -K \sqrt{\\frac{4D}{\pi}} \sqrt{t}$', ha='left', va='bottom', color='red', fontsize=20)
ax[4].text(10, 0, '$Δc = -K \sqrt{\\frac{4D t_i}{\pi}} \\left[ \\frac{\sqrt{t_i}}{\sqrt{t} + \sqrt{t - t_i} } \\right]$', ha='left', va='bottom', color='red', fontsize=20)

# INSETS
x = np.linspace(0, 10, 200)
ax3_in = ax[3].inset_axes([0.5, 0.05, 0.48, 0.48])
t_ = np.logspace(-2, 1)
delta_C = func_x0(K, D, t_)
ax3_in.plot(np.sqrt(t_), delta_C, 'ko')

ax3_in.xaxis.tick_top()
ax3_in.set_xlim(xmin=0)
ax3_in.set_ylim(ymin=0)
ax3_in.set_ylabel('$\mathregular{Δc \ (x = 0)}$', fontsize=15)
ax3_in.set_xlabel('$\mathregular{\sqrt{t}}$', fontsize=15)
ax3_in.xaxis.set_label_position('top')





ax4_in = ax[4].inset_axes([0.37, 0.37, 0.58, 0.58])
t_ = np.logspace(1, 2)
ti = 10
tau = np.sqrt(ti)/(np.sqrt(t_) + np.sqrt(t_ - ti))
delta_C = func_x0(K, D, t_, ti=ti)

ax4_in.plot(tau, delta_C, 'ko')

ax4_in.set_xlim(xmin=-0.02)
ax4_in.set_ylim(ymin=-0.1)
ax4_in.set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=15)
ax4_in.set_xlabel('$\mathregular{\\tau}$', fontsize=15)
ax4_in.text(1.0, 0.0, '$\mathregular{\\tau = \\frac{\sqrt{t_i}}{\sqrt{t} + \sqrt{t - t_i} } }$', ha = 'right', va = 'bottom', fontsize=20, color = 'red')

fig.tight_layout()

os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Solid-State-Electrolyte/Model/Equations/')
fig.savefig('MathStackExchangeAnswer.jpg', format='jpg', dpi=300, bbox_inches = "tight")

# In[2]:
'''
    When the current is on, the expression for the concentration can be modified, and a transformed variable can be used
    c(x,t) = c_∞ - K[ √(4Dt/π) exp(-x^2/(4Dt)) - x erfc(x/(√(4Dt))) ]
    c(x,t) - c_∞ = Δc
    -1/K ⋅ Δc/(√(4Dt)) = 1/√π exp(-x^2/(4Dt)) - x/√(4Dt) erfc(x/(√(4Dt)))
    -1/K ⋅ Δc/(√(4Dt)) = 1/√π exp(-η^2) - η erfc(η)     ; η = x/√(4Dt)
'''

rows, columns = [1, 1]
ax, fig = axes(1, rows, columns)

eta = np.logspace(-5,np.log10(2.1),200)
mod_delC = np.sqrt(1/np.pi) * np.exp(-eta**2) - eta * scipy.special.erfc(eta)

x = 1.5
print( (np.sqrt(1/np.pi) * np.exp(-x**2) - x * scipy.special.erfc(x)) * 100)

ax[1].plot(eta, mod_delC)
ax[1].set_xlabel('$\mathregular{η = \\frac{x}{\sqrt{4Dt}}}$', fontsize=20)
# ax[1].set_ylabel('$\mathregular{\\frac{Δc}{K\sqrt{4Dt}}}$', fontsize=20)
ax[1].set_ylabel('$\mathregular{\\frac{z_+ \\nu_+ F D}{-(1 - t_+^0)I} \\frac{Δc}{\sqrt{4Dt}}}$', fontsize=20)



# In[2]:
# Import 
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/')

data_TP_model = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])
print(data_TP_model.keys())


# In[3]:
'''
    Comparison between
        1. analytical solution whihch assumes an infinite diffusion length, i.e. c(t,x=∞) = c_bulk
            (solid colored lines)
        2. Numerical simulation where length is not ∞ (x = 2.0 cm)

Conclusion:
    Despite a slight difference in assumptions, the results are nearly identical
        validation of the infinite diffusion length assumption
'''

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)


D = 1e-5
trans_Li = 0.5
z_Li = +1
nu_Li = 1
i_applied_cm2 = 10e-6
K = -(1 - trans_Li)*i_applied_cm2/(z_Li * nu_Li * Fconst * D)
x_hat = np.logspace(-5, np.log10(0.3),100)
x_hat = np.insert(x_hat, 0, 0) # add zero to beginning of array
ti = 180

for time in [0, 0.1, 10, 50, ti]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[1].plot(data_['Position'], data_['Delta_Conc'], linewidth = 4)
    ax[1].plot(data_['Position'], func(K, D, time, data_['Position']), '--k', linewidth = 1.5)
    ax[1].plot(x_hat, func(K, D, time, x_hat), '--k', linewidth = 1.5)


for time in [ti, 190, 230, 360, data_TP_model['Time'].max()]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[2].plot(data_['Position'], data_['Delta_Conc'], linewidth = 4)
    ax[2].plot(x_hat, func(K, D, time, x_hat, ti=ti), '--k', linewidth = 1.5)

ax[1].set_xlim([-1e-2, 0.2])
ax[2].set_xlim([-1e-2, 0.2])



# x = 0 data
data_ = data_TP_model[data_TP_model['Position'] == 0.0]
data_1 = data_[data_['Time'] <= 180]
data_2 = data_[data_['Time'] >= 180]
ax[3].plot(data_1['Time'], data_1['Delta_Conc'], linewidth = 4, color='darkorange')
ax[4].plot(data_2['Time'], data_2['Delta_Conc'], linewidth = 4, color='darkorange')

ylim = ax[3].get_ylim()
ax[4].set_ylim(ylim)


ax[3].plot(data_1['Time'], func_x0(K, D, data_1['Time']), '--k', linewidth = 1.5)
ax[4].plot(data_2['Time'], func_x0(K, D, data_2['Time'], ti=ti), '--k', linewidth = 1.5)

ax[1].set_ylabel('$\mathregular{\Delta c}$', fontsize=20)
ax[3].set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=20)
ax[2].set_yticklabels([])
ax[4].set_yticklabels([])

ax[1].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[2].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[3].set_xlabel('$t$', fontsize=20)
ax[4].set_xlabel('$t$', fontsize=20)

fig.tight_layout()

fig.savefig('TransferencePolarization_AnalyticalSoln_vs_FiniteSimultion.pdf', format='pdf', dpi=300, bbox_inches = "tight")

# In[5]:
# Dilute Solution Theory, symmetric cell for transference polarization
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/Symmetric/')

data_TP_model_sym = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])
print(data_TP_model_sym.keys())

D_Li = 2e-6
D_PF6 = 1e-6
t_Li = (D_Li)/(D_Li + D_PF6)
z_Li = +1.0
z_PF6 = -1.0
D_eff = 2*D_Li*D_PF6/(z_Li*D_Li - z_PF6*D_PF6)
# D_eff = (z_+ u_+ D_- - z_- u_- D_+) / (z_+ u_+ - z_- u_-)
# u_i = D_i / (RT)
print(D_eff)
c_bulk = 1e-3 # mol/cm3
L = 1.0 # cm

print("Max Diffusion Rate (mA/cm2)", 2*D_eff*c_bulk/L*Fconst * 1e3)

# F∇Φ = -(D_+ - D_-)/(z_+ u_+ - z_- u_-) ∇lnc = F∇Φ = -(D_+ - D_-)/(z_+ u_+ - z_- u_-) 1/c ∇c
# ΔΦ = -(D_+ - D_-)/(z_+ u_+ - z_- u_-) Δln(c) / F

# In[6]:
# Voltage Profiles

ax, fig = axes(1, rows=2, columns=3)

edge_ = data_TP_model_sym['Position'] == max(data_TP_model_sym['Position'])
data_TP_model_sym_edge = data_TP_model_sym[edge_]

new_index = range(len(data_TP_model_sym_edge))
data_TP_model_sym_edge = data_TP_model_sym_edge.reset_index(drop=True)

Pot_NJ = data_TP_model_sym_edge['Potential']
time = data_TP_model_sym_edge['Time'] / 3600.

index = data_TP_model_sym_edge.index
relaxation = data_TP_model_sym_edge['Current'] == 0.0
relax_indices = index[relaxation]
relax_indices = relax_indices.tolist()

polarization = data_TP_model_sym_edge['Current'] != 0.0
polar_indices = index[polarization]
polar_indices = polar_indices.tolist()

rest_list = []
for k, g in groupby(enumerate(relax_indices), lambda ix : ix[0] - ix[1]):
    rest_list.append(list(map(itemgetter(1), g)))

polar_list = []
for k, g in groupby(enumerate(polar_indices), lambda ix : ix[0] - ix[1]):
    polar_list.append(list(map(itemgetter(1), g)))    

# #
ax[1].plot(time, abs(Pot_NJ)*1e3, 'b-', markersize=1)
ax[1].plot(time[relaxation], abs(Pot_NJ)[relaxation]*1e3, 'r.', markersize=1)

ax[2].plot(time[relaxation], abs(Pot_NJ)[relaxation]*1e3, 'r.', markersize=1)

Conc_0  = data_TP_model_sym[data_TP_model_sym['Position'] == 0]['Conc'].values
Conc_NJ = data_TP_model_sym[data_TP_model_sym['Position'] == max(data_TP_model_sym['Position'])]['Conc'].values
# ΔΦ = -(D_+ - D_-)/(z_+ u_+ - z_- u_-) Δln(c) / F
# ΔΦ = -(D_+ - D_-)/(z_+ D_+ - z_- D_-) Δln(c) * (RT/F)
Delta_Phi = -(D_Li - D_PF6)/(z_Li * D_Li - z_PF6*D_PF6) * Rigc*Temp/Fconst * (np.log(Conc_0) - np.log(Conc_NJ))
ax[2].plot(time, abs(Delta_Phi)*1e3, 'k-', zorder=-100)


conc_hat = np.logspace(-3, np.log10(2)) # mol/L (M)
Phi_hat = -(D_Li - D_PF6)/(z_Li * D_Li - z_PF6*D_PF6) * Rigc*Temp/Fconst * np.log(conc_hat/1e3)
ax[3].plot(conc_hat, Phi_hat*1e3)
ax[3].set_xlabel('Concentration (M)')
ax[3].set_ylabel('U (mV)')
#
I_sqrt_ti = []
del_U_list = []

for o, lst in enumerate(rest_list):
    delta_t = time[lst] - time[lst[0]-1]
    if o == 0:
        t_i = time[lst[0]] - 0.0
    else:
        t_i = time[lst[0]] - time[rest_list[o-1][-1]]
    t = delta_t + t_i
    tau = np.sqrt(t_i)/(np.sqrt(t) + np.sqrt(t - t_i))
    Pot = Pot_NJ[lst]

    I_app = data_TP_model_sym_edge['Current'][lst[0]-1]

    ax[4].plot(delta_t*60, Pot*1e3, '-')

    ax[5].plot(tau, Pot*1e3, '-')

    # fit_window = (tau < 0.5)
    fit_window = (tau > 0.3) & (tau < 0.6)
    p = np.polyfit(tau[fit_window], Pot[fit_window]*1e3, 1)

    tau_hat = np.linspace(0,1)
    ax[5].plot(tau_hat, np.polyval(p, tau_hat), 'k--', linewidth = 1.5)

    del_U = p[0] #np.polyval(p, 1) - np.polyval(p, 0) - slope
    ax[6].plot(I_app*np.sqrt(t_i)*1e3, del_U, 'ks')

    I_sqrt_ti.append(I_app*np.sqrt(t_i)*1e3)
    del_U_list.append(del_U)

#
ax[1].set_ylabel('Measured Potential (mV)')
ax[1].set_xlabel('Time (hours)')

ax[2].set_ylabel('Measured Potential (mV)')
ax[2].set_xlabel('Time (hours)')

ax[4].set_ylabel('Measured Potential (mV)')
ax[4].set_xlabel('Δt (minutes)')

ax[5].set_ylabel('Measured Potential (mV)')
ax[5].set_xlabel('τ')
ax[5].set_xlim(xmin=0)
ax[5].set_ylim(ymin=0)

ax[6].set_ylabel('ΔU (mV)')
ax[6].set_xlabel('$\mathregular{I t_i^{1/2} \ (mA \ cm^{-2} \ s^{1/2})}$')
ax[6].set_xlim(xmin=0)
ax[6].set_ylim(ymin=0)

p = np.polyfit(I_sqrt_ti, del_U_list, 1)
I_hat = np.linspace(0, max(I_sqrt_ti))
ax[6].plot(I_hat, np.polyval(p, I_hat), 'k--', linewidth=1)

fig.tight_layout()
fig.savefig('TransferPolar_SymmDiluteSolnMigr_vs_AnalyticSoln__Potentials.pdf', format='pdf', dpi=300, bbox_inches = "tight")
# fig.savefig('TransferPolar_SymmDiluteSolnMigr_vs_AnalyticSoln_UniformConc_Potentials.pdf', format='pdf', dpi=300, bbox_inches = "tight")

# In[7]:
# Concentration Profiles

rows, columns = [2, 3]
ax, fig = axes(1, rows, columns)

initial_conc = data_TP_model_sym['Conc'].values[0]
delta_C = data_TP_model_sym['Conc'].values - initial_conc
data_TP_model_sym['Delta_Conc'] = delta_C

data_ = data_TP_model_sym[data_TP_model_sym['Position'] == 0.0]
ax[1].plot(data_['Time']/3600., data_['Delta_Conc'])

for o, lst in enumerate(polar_list):
    delta_time = data_['Time'].values[lst] - data_['Time'].values[lst[0]]

    I_app = data_['Current'].values[lst[0]]
    K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)

    ax[2].plot(delta_time, data_['Delta_Conc'].values[lst])
    ax[2].plot(delta_time, func_x0(K, D_eff, delta_time, ti=np.inf), 'k--', linewidth=2)

    ax[5].plot(np.sqrt(delta_time), data_['Delta_Conc'].values[lst])
    ax[5].plot(np.sqrt(delta_time), func_x0(K, D_eff, delta_time, ti=np.inf), 'k--', linewidth=2)

x_inset_pos, y_inset_pos = (0.35, 0.35)
pad = 0.04
ax3_in = ax[3].inset_axes([x_inset_pos, y_inset_pos, 1. - x_inset_pos - pad, 1. - y_inset_pos - pad])
# x_inset_pos, y_inset_pos = (0.35, 0.35)
x_pos = 0.12
y_pos = 0.65
ax6_in = ax[6].inset_axes([x_pos, y_pos, 0.57-x_pos, 1 - y_pos - 0.02])
# ax6_in.xaxis.tick_top()
# ax6_in.xaxis.set_label_position('top')
ax6_in.set_xlabel('τ')
ax6_in.set_ylabel('Φ (mV)')
ax6_in.set_xlim((-0.03, 1.03 ))


for o, lst in enumerate(rest_list):
    delta_time = data_['Time'].values[lst] - data_['Time'].values[polar_list[o][-1]]
    print(delta_time[0])
    delta_time_hat = np.logspace(-3,np.log10(max(delta_time)))
    delta_time_hat = np.insert(delta_time_hat, 0, 0)
    
    I_app = data_TP_model_sym_edge['Current'][lst[0]-1]
    K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)
    ti = 60

    ax[3].plot(delta_time, data_['Delta_Conc'].values[lst])
    ax[3].plot(delta_time_hat, func_x0(K, D_eff, delta_time_hat + ti, ti=ti), 'k--', linewidth=2)

    ax3_in.plot(delta_time, data_['Delta_Conc'].values[lst])
    ax3_in.plot(delta_time_hat, func_x0(K, D_eff, delta_time_hat + ti, ti=ti), 'k--', linewidth=2)    

    t = delta_time + ti
    tau = np.sqrt(ti)/(np.sqrt(t) + np.sqrt(t - ti))
    tau_hat = np.linspace(0,1)

    ax[6].plot(tau, data_['Delta_Conc'].values[lst])
    ax[6].plot(tau_hat, -K*np.sqrt(4*D_eff*ti/np.pi)*tau_hat, 'k--', linewidth=2)
    ax6_in.plot(tau, Pot_NJ[lst]*1e3)


    # get slope
    fit_window = (tau > 0.3) & (tau < 0.6)
    # fit_window = (tau < 0.6)
    p = np.polyfit(tau[fit_window], data_['Delta_Conc'].values[lst][fit_window], 1)
    ax[6].plot(tau_hat, np.polyval(p, tau_hat), 'r-', linewidth=2)
    # print(p[0], -K*np.sqrt(4*D_eff*ti/np.pi))
    ax[4].plot(I_app*np.sqrt(ti)*1e3, -K*np.sqrt(4*D_eff*ti/np.pi), 'ks')
    ax[4].plot(I_app*np.sqrt(ti)*1e3, p[0], 'ro')

ax3_in.set_xlim(-1, 25)
ax3_in.set_xlabel('Time (seconds)')
ax3_in.set_ylabel('Δc (mol/cm3)')

for i in range(1, rows*columns+1):
    ax[i].set_ylabel('Δc (mol/cm3)')

ax[1].set_xlabel('Time (hours)')
ax[2].set_xlabel('Time (seconds)')
ax[3].set_xlabel('Time (seconds)')

ax[4].set_xlabel('$\mathregular{I t_i^{1/2} \ (mA \ cm^{-2} \ s^{1/2})}$')
ax[5].set_xlabel('$\mathregular{\sqrt{t} \ (s^{1/2})}$')
ax[6].set_xlabel('τ (-)')

fig.tight_layout()
fig.savefig('TransferPolar_SymmDiluteSolnMigr_vs_AnalyticSoln.pdf', format='pdf', dpi=300, bbox_inches = "tight")
# fig.savefig('TransferPolar_SymmDiluteSolnMigr_vs_AnalyticSoln_UniformConc.pdf', format='pdf', dpi=300, bbox_inches = "tight")

# In[8]:
'''
    Import the LiPF6 in EC:EMC Concentrated Solution Theory Model Data 
    for Transference Polarization
'''
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_LiPF6_in_EC_EMC/PseodoExpts/TransferPolar/')

data_files = glob.glob("Time_Voltage_Position*.txt")
data_files.sort(key=natural_keys)
print(data_files)

# data_TP_conc_soln_model = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])

# In[9]:
print(data_TP_conc_soln_model.keys())
# t_Li, D_eff = ( 0.38608705750000089,        5.3471515528662765e-006) # c_bulk = 1e-5
t_Li, D_eff = ( 0.37142806000000156,        5.0550283571730083E-006) # c_bulk = 1e-4
t_Li, D_eff = ( 0.30332030000000043,        3.9380756371031145E-006) # c_bulk = 5e-4
t_Li, D_eff = ( 0.21919360000000010,        2.8822697547954958E-006) # c_bulk = 1e-3
t_Li, D_eff = ( 0.10417520000000291,        1.5439590128016643E-006) # c_bulk = 2e-3

t_Li_D_eff_list = (
    ( 0.37142806000000156,        5.0550283571730083E-006),
    ( 0.30332030000000043,        3.9380756371031145E-006),
    ( 0.21919360000000010,        2.8822697547954958E-006),
    ( 0.10417520000000291,        1.5439590128016643E-006),
)

K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)

for o, data_file in enumerate(data_files):

    t_Li, D_eff = t_Li_D_eff_list[o]
    print(data_file, t_Li, D_eff, o)

    file_name, ext = os.path.splitext(data_file)
    conc_text = file_name[-4:]
    print(file_name, conc_text)


# In[10]:
# plot concentrated solution data

for o, data_file in enumerate(data_files[:]):
    # Import the data file as a dataframe
    data_TP_conc_soln_model = pd.read_table(data_file, delim_whitespace=True, header=[0], skiprows=[1])

    # the bulk concentration transport properties
    t_Li, D_eff = t_Li_D_eff_list[o]
    
    file_name, ext = os.path.splitext(data_file)
    conc_text = file_name[-4:]  # the bulk concentration value in mol/L

    print(data_file, t_Li, D_eff, conc_text)

    # get initial bulk concentration values
    initial_conc = data_TP_conc_soln_model['c_LiPF6'].values[0]
    # dev_C = c_edge - c_bulk
    dev_C = data_TP_conc_soln_model['c_LiPF6'].values - initial_conc
    data_TP_conc_soln_model['Dev_Conc'] = dev_C

    # the values of concentration and potential at the system boundaries
    data_edges = data_TP_conc_soln_model[data_TP_conc_soln_model['Position'] == 0.0]
    data_NJ = data_TP_conc_soln_model[data_TP_conc_soln_model['Position'] == max(data_TP_conc_soln_model['Position'])]

    # The indices for current OFF
    new_index = range(len(data_edges))
    data_edges = data_edges.reset_index(drop=True)
    data_edges = data_edges.rename(columns={'c_LiPF6': 'c_LiPF6_x=0'})
    data_edges['c_LiPF6_x=NJ'] = data_NJ['c_LiPF6'].values
    data_edges['Delta_C'] = data_edges['c_LiPF6_x=0'] - data_edges['c_LiPF6_x=NJ']

    index = data_edges.index
    relaxation = data_edges['Current'] == 0.0
    relax_indices = index[relaxation]
    relax_indices = relax_indices.tolist()

    rest_list = []
    for k, g in groupby(enumerate(relax_indices), lambda ix : ix[0] - ix[1]):
        rest_list.append(list(map(itemgetter(1), g)))

    # indices of current ON
    polarization = data_edges['Current'] != 0.0
    polar_indices = index[polarization]
    polar_indices = polar_indices.tolist()

    polar_list = []
    for k, g in groupby(enumerate(polar_indices), lambda ix : ix[0] - ix[1]):
        polar_list.append(list(map(itemgetter(1), g)))


    # Plot
    # 1) c_edge - c_bulk vs experiment time     2) c_edge - c_bulk vs time (current ON)     3) c_edge - c_bulk vs time (current OFF)
    # 4) c_edge - c_bulk vs I√tᵢ                5) c_edge - c_bulk vs √t                    4) c_edge - c_bulk vs τ
    rows, columns = [2, 3]
    ax, fig = axes(1+o, rows, columns)

    data_ = data_edges
    # conc_data = data_['Dev_Conc']
    conc_data = data_['Delta_C']/2
    ax[1].plot(data_['Time'], conc_data, 'k-')

    for o, lst in enumerate(polar_list):
        delta_time = data_['Time'].values[lst] - data_['Time'].values[lst[0]]
        delta_time = delta_time*3600.
        delta_time_hat = np.logspace(-3,np.log10(max(delta_time)))
        delta_time_hat = np.insert(delta_time_hat, 0, 0)

        I_app = data_['Current'].values[lst[0]]
        K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)

        line, = ax[2].plot(delta_time, conc_data.values[lst])
        ax[2].plot(delta_time_hat, func_x0(K, D_eff, delta_time_hat, ti=np.inf), 'k--', linewidth=2)
        color = line.get_color() 

        ax[1].plot(data_['Time'].values[lst], conc_data.values[lst], '.', color=color)

        ax[5].plot(np.sqrt(delta_time), conc_data.values[lst])
        ax[5].plot(np.sqrt(delta_time_hat), func_x0(K, D_eff, delta_time_hat, ti=np.inf), 'k--', linewidth=2)


    x_inset_pos, y_inset_pos = (0.35, 0.35)
    pad = 0.08
    ax3_in = ax[3].inset_axes([x_inset_pos, y_inset_pos, 1. - x_inset_pos - pad, 1. - y_inset_pos - pad])

    for o, lst in enumerate(rest_list):
        delta_time = data_['Time'].values[lst] - data_['Time'].values[lst[0]]
        delta_time = delta_time*3600.
        delta_time_hat = np.logspace(-3,np.log10(max(delta_time)))
        delta_time_hat = np.insert(delta_time_hat, 0, 0)
        
        # print(data_)
        # print(len(lst), lst[0])

        I_app = data_['Current'].values[lst[0]-1]
        K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)
        ti = 60

        line, = ax[3].plot(delta_time, conc_data.values[lst])
        ax[3].plot(delta_time_hat, func_x0(K, D_eff, delta_time_hat + ti, ti=ti), 'k--', linewidth=2, zorder=100)
        ax3_in.plot(delta_time, conc_data.values[lst])
        ax3_in.plot(delta_time_hat, func_x0(K, D_eff, delta_time_hat + ti, ti=ti), 'k--', linewidth=2, zorder=100)

        color = line.get_color() 
        ax[1].plot(data_['Time'].values[lst], conc_data.values[lst], '.', color=color)

        t = delta_time + ti
        tau = np.sqrt(ti)/(np.sqrt(t) + np.sqrt(t - ti))
        tau_hat = np.linspace(0,1)

        ax[6].plot(tau, conc_data.values[lst])
        ax[6].plot(tau_hat, -K*np.sqrt(4*D_eff*ti/np.pi)*tau_hat, 'k--', linewidth=2, zorder=200)

        # get slope
        fit_window = (tau > 0.3) &(tau < 0.7)
        p = np.polyfit(tau[fit_window], conc_data.values[lst][fit_window], 1)
        ax[6].plot(tau_hat, np.polyval(p, tau_hat), 'r-', linewidth=2, zorder=100)
        # print(p[0], -K*np.sqrt(4*D_eff*ti/np.pi))
        
        ax[4].plot(I_app*np.sqrt(ti)*1e3, p[0], 'ro')

        K = -(1-t_Li)*I_app/(z_Li*Fconst*D_eff)
        ax[4].plot(I_app*np.sqrt(ti)*1e3, -K*np.sqrt(4*D_eff*ti/np.pi), 'ks')
        # print(p[0], -K*np.sqrt(4*D_eff*ti/np.pi), p[0] / (-K*np.sqrt(4*D_eff*ti/np.pi)) )

    for i in range(1, rows*columns+1):
        ax[i].set_ylabel('$\mathregular{\Delta c \ (mol/cm^3)}$')

    ax[1].set_title(conc_text + ' $\mathregular{LiPF_6}$ in EC:EMC', fontsize=20, fontweight='bold')
    ax[2].set_title('Current ON Data',  fontsize=20, fontweight='bold')
    ax[3].set_title('Current OFF Data', fontsize=20, fontweight='bold')

    ax[1].set_xlabel('Time (hours)')
    ax[2].set_xlabel('Time (seconds)')
    ax[3].set_xlabel('Time (seconds)')
    ax[3].set_xlim(-50, 1800)
    ax3_in.set_xlim(-1, 25)
    ax3_in.set_xlabel('Time (seconds)')
    ax3_in.set_ylabel('$\mathregular{\Delta c \ (mol/cm^3)}$')

    ax[4].set_xlabel('$\mathregular{I t_i^{1/2} \ (mA \ cm^{-2} \ s^{1/2})}$')
    ax[5].set_xlabel('$\mathregular{\sqrt{t} \ (s^{1/2})}$')
    ax[6].set_xlabel('τ (-)')    

    fig.tight_layout()
    fig.savefig('TransferPolar_SymmConcSoln_vs_AnalyticSoln_{}.pdf'.format(conc_text), 
                format='pdf', dpi=300, bbox_inches = "tight")
    # fig.savefig('TransferPolar_SymmConcSoln_vs_AnalyticSoln_UniformConc_{}.pdf'.format(conc_text), 
    #             format='pdf', dpi=300, bbox_inches = "tight")
# In[12]:
# Concentration Cell
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/26_LiPF6_in_EC_EMC/PseodoExpts/ConcentrationCell/')
data_conc_soln_ConcCell = pd.read_table('Time_Voltage_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])
# In[13]:

unique_time = np.unique(data_conc_soln_ConcCell['Time'].values)

data_edges = data_conc_soln_ConcCell[data_conc_soln_ConcCell['Position'] == 0.0]
data_NJ = data_conc_soln_ConcCell[data_conc_soln_ConcCell['Position'] == max(data_conc_soln_ConcCell['Position'])]

rows, columns = [1,1]
ax, fig = axes(1, rows, columns)
conc_data = data_conc_soln_ConcCell[data_conc_soln_ConcCell['Time'] == unique_time[10]]['c_LiPF6']*1e3
potential = data_conc_soln_ConcCell[data_conc_soln_ConcCell['Time'] == unique_time[10]]['Φ_2']*1e3
ax[1].semilogx(conc_data, potential - max(potential), 'ko-')

ax[1].set_xlim(0.08, 4)
ax[1].set_ylim(-300, 10)
ax[1].set_ylabel('$\mathregular{\Phi_2 \ (mV)}$')
ax[1].set_xlabel('$\mathregular{c_{LiPF_6} \ (mol/L)}$')

fig.tight_layout()

# In[100]:
'''
    Restricted Diffusion Experiments
'''

def time_window_return_P(data_time, data_conc, low_limit, upper_limit):
    moving_avg_window = 60
    delta_time_slice = data_time[(data_time >= low_limit) & (data_time <= upper_limit)]
    voltage_slice = data_conc[(data_time >= low_limit) & (data_time <= upper_limit)]

    p = np.polyfit(delta_time_slice[moving_avg_window-1:],
                    np.log(voltage_slice.rolling(moving_avg_window).mean().dropna()), 1)

    return p

for percent in [1, 0.5, 0.1, 0.01]:
    os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/19_RestrictedDiffusion/Normal/1/')

    data = read_Fortran('Time_Conc_Position.txt')

    rows, columns = [4,2]
    ax, fig = axes(1, rows, columns)

    data_edge = data[data['Position'] == 0]
    ax[1].plot(data_edge['Time'], -data_edge['Delta_Conc']*1e3, '-', color='black', zorder=100)
    ax[2].plot(data_edge['Time'], np.log(-data_edge['Delta_Conc']*1e3), '-', color='black', zorder=100)

    ''' With Random Noise '''
    random_noise = (np.random.rand(len(data_edge))*2 - 1) * 2e-1
    real_data = -data_edge['Delta_Conc']*1e3 + random_noise
    data_edge['Real_Data'] = real_data.values
    ax[1].plot(data_edge['Time'], real_data, 'o', color='lightblue', zorder=-10)
    ax[2].plot(data_edge['Time'], np.log(real_data), 'o', color='lightblue', zorder=-10)

    ''' Moving Average '''
    ax[1].plot(data_edge['Time'], real_data.rolling(60).mean(), 'o', color='blue', zorder=2)
    ax[2].plot(data_edge['Time'], np.log(real_data.rolling(60).mean()), 'o', color='blue', zorder=-10)

    ax[1].set_xlim([-1, 20])
    ax[2].set_xlim([-1, 10])
    ax[2].set_ylim([-5, 0.5])

    ''' Get the Diffusion Coefficient '''
    ''' Exact Data '''
    data_fit = data_edge[data_edge['Time'] <= 20]
    data_fit = data_fit[data_fit['Time'] >= 1]
    p = np.polyfit(data_fit['Time'], np.log(-data_fit['Delta_Conc']*1e3), 1)
    print("Exact:", -p[0] / 3600 * 1**2 / np.pi**2)

    ''' Real Noisy Data '''
    p = time_window_return_P(data_edge['Time'], data_edge['Real_Data'], 1, 10)
    time_hat = np.linspace(0,20)
    ax[2].plot(time_hat, np.polyval(p, time_hat), '-', color='darkblue')
    print("Real:", -p[0] / 3600 * 1**2 / np.pi**2)

    fig.suptitle('Restricted Diffusion', y = 0.93)

    ax[1].set_ylabel('Δc  (M)')
    ax[1].set_xlabel('Time (hours)')
    ax[2].set_ylabel('ln(Δc)')
    ax[2].set_xlabel('Time (hours)')

    fig.tight_layout()
#


# In[101]:
rows, columns = [4,2]
ax, fig = axes(1, rows, columns)

for o, percent in enumerate([1, 0.5, 0.1, 0.01]):
    os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/19_RestrictedDiffusion/Normal/')
    os.chdir(str(percent))

    data = read_Fortran('Time_Conc_Position.txt')

    data_edge = data[data['Position'] == 0]
    ax[1 + o*2].plot(data_edge['Time'], -data_edge['Delta_Conc']*1e3, '-', color='black', zorder=100)
    ax[2 + o*2].plot(data_edge['Time'], np.log10(-data_edge['Delta_Conc']*1e3), '-', color='black', zorder=100)

    ''' With Random Noise '''
    random_noise = (np.random.rand(len(data_edge))*2 - 1) * 1e-2
    real_data = -data_edge['Delta_Conc']*1e3 + random_noise
    data_edge['Real_Data'] = real_data.values
    ax[1 + o*2].plot(data_edge['Time'], real_data, 'o', color='lightblue', zorder=-10)
    ax[2 + o*2].plot(data_edge['Time'], np.log10(real_data), 'o', color='lightblue', zorder=-10)

    ''' Moving Average '''
    ax[1 + o*2].plot(data_edge['Time'], real_data.rolling(60).mean(), 'o', color='blue', zorder=2)
    ax[2 + o*2].plot(data_edge['Time'], np.log10(real_data.rolling(60).mean()), 'o', color='blue', zorder=-10)

    ax[1 + o*2].set_xlim([-1, 20])
    ax[2 + o*2].set_xlim([-1, 10])
    if o == 0:
        ax[2 + o*2].set_ylim([-2, 0.])
    elif o == 1:
        ax[2 + o*2].set_ylim([-2, 0.])
    elif o == 2:
        ax[2 + o*2].set_ylim([-3, -1])
    elif o == 3:
        ax[2 + o*2].set_ylim([-4, -1.5])

    ''' Get the Diffusion Coefficient '''
    ''' Exact Data '''
    data_fit = data_edge[data_edge['Time'] <= 20]
    data_fit = data_fit[data_fit['Time'] >= 1]
    p = np.polyfit(data_fit['Time'], np.log(-data_fit['Delta_Conc']*1e3), 1)
    print("Exact:", -p[0] / 3600 * 1**2 / np.pi**2)

    ''' Real Noisy Data '''
    p = time_window_return_P(data_edge['Time'], data_edge['Real_Data'], 1, 6)
    time_hat = np.linspace(0,20)
    ax[2 + o*2].plot(time_hat, np.polyval(p/np.log(10), time_hat), '-', color='darkblue')
    print("Real:", -p[0] / 3600 * 1**2 / np.pi**2)

    ax[1 + o*2].set_ylabel('Δc  (M)')
    ax[1 + o*2].set_xlabel('Time (hours)')
    ax[2 + o*2].set_ylabel('$\mathregular{log_{10}}$(Δc)')
    ax[2 + o*2].set_xlabel('Time (hours)')

fig.suptitle('Restricted Diffusion')
fig.tight_layout()

# In[10]:
'''
    Situation where solvent is not charge neutral and contributes to the ionic current
'''

def Nernst_Potential(c_Li, c_sol):
    # c_Li - lithium concentration
    # c_sol - solvent concentration

    Rigc = 8.314
    Temp = 300
    Fconst = 96485
    U_ref = 0

    return U_ref - Rigc * Temp / Fconst * np.log(c_Li / c_sol)

def Nernst_Potential_theta_max(c_Li, c_sol):
    # c_Li - lithium concentration
    # c_sol - solvent concentration

    Rigc = 8.314
    Temp = 300
    Fconst = 96485
    U_ref = 0

    theta = c_Li / (4*c_sol)

    return U_ref - Rigc * Temp / Fconst * np.log(theta/(1-theta) )

rows, columns = [1,1]
ax, fig = axes(1, rows, columns)

c_Li = np.linspace(1e-3, 4-1e-3, 100)
c_sol = np.ones(len(c_Li))*2

ax[1].plot(c_Li, Nernst_Potential(c_Li, c_sol))
ax[1].plot(c_Li, Nernst_Potential_theta_max(c_Li, c_sol))

ax[1].set_ylim(-0.1, 0.1)

# In[11]:

t_Li = 0.1
t_EMI = 0.5

K_Li = -1
K_EMI = K_Li * - t_EMI / (1 - t_Li)

t = np.linspace(0, 2)

rows, columns = [1,2]
ax, fig = axes(1, rows, columns)

delta_C_Li = func_x0(K_Li, D, t)
delta_C_EMI = func_x0(K_EMI, D, t)

ax[1].plot(t, delta_C_Li)
ax[1].plot(t, delta_C_EMI)