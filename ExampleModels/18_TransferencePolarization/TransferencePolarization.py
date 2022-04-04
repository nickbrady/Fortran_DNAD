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

def read_Fortran(file):
    return pd.read_table(file, delim_whitespace=True, header=[0], skiprows=[1])


# def axes(number, rows, columns):
#     fig = plt.figure(number, figsize=(6*columns, 5*rows), dpi=80)
#     ax = []
#     ax.append([])
#     for i in range(1, rows*columns+1):
#         ax.append(fig.add_subplot(rows, columns, i))
#
#     return ax, fig

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

rows, columns = [2, 2]
ax, fig = axes(1, rows=rows, columns=columns)

def func(K, D, t, x, ti=False):
    def f(D, t, x):
        return np.sqrt(4*D*t/np.pi) * np.exp(-x**2/(4*D*t)) - x * scipy.special.erfc(x/np.sqrt(4*D*t))


    f_less_ti = f(D, t, x)

    f_gt_ti = 0.0


    if ti:
        if (t > ti):
            f_gt_ti = f(D, t - ti, x)

    return K * (f_less_ti - f_gt_ti)


def func_x0(K, D, t, ti=False):

    if not(ti):
        ti = np.inf


    delta_c = np.where(t < ti, np.sqrt(4*D*t/np.pi), np.sqrt(4*D*t/np.pi) - np.sqrt(4*D*(t-ti)/np.pi))

    return delta_c*K


K = 1
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

    ax[1].text(10, 3.5, '$t < t_i$', va='top', ha='right', fontsize=20)
    ax[2].text(10, 3.5, '$t > t_i$', va='top', ha='right', fontsize=20)


    ax[3].set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=20)

    ax[3].set_xlabel('$t$', fontsize=20)
    ax[4].set_xlabel('$t$', fontsize=20)

    ax[3].text(0,  3.5, '$t < t_i$', va='top', ha='left', fontsize=20)
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

ax[3].text(0, 2.5, '$\Delta c = K \sqrt{\\frac{4D}{\pi}} \sqrt{t}$', ha='left', va='bottom', color='red', fontsize=20)
ax[4].text(10, 0, '$\Delta c = K \sqrt{\\frac{4D t_i}{\pi}} \\left[ \\frac{\sqrt{t_i}}{\sqrt{t} + \sqrt{t - t_i} } \\right]$', ha='left', va='bottom', color='red', fontsize=20)

# INSET
K = 1
x = np.linspace(0, 10, 200)
ax3_in = ax[3].inset_axes([0.5, 0.05, 0.48, 0.48])
t_ = np.logspace(-2, 1)
for t in t_:
    delta_C = func(K, D, t, x)
    ax3_in.plot(np.sqrt(t), delta_C[0], 'ko')

ax3_in.xaxis.tick_top()
ax3_in.set_xlim(xmin=0)
ax3_in.set_ylim(ymin=0)
ax3_in.set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=15)
ax3_in.set_xlabel('$\mathregular{\sqrt{t}}$', fontsize=15)
ax3_in.xaxis.set_label_position('top')





ax4_in = ax[4].inset_axes([0.37, 0.37, 0.58, 0.58])
t_ = np.logspace(1, 2)
ti = 10
tau = np.sqrt(ti)/(np.sqrt(t_) + np.sqrt(t_ - ti))
for t in t_:
    delta_C = func(K, D, t, x, ti=ti)
    tau = np.sqrt(ti)/(np.sqrt(t) + np.sqrt(t - ti))
    ax4_in.plot(tau, delta_C[0], 'ko')

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
    c(x,t) = c_∞ + K[ √(4Dt/π) exp(-x^2/(4Dt)) - x erfc(x/(√(4Dt))) ]
    c(x,t) - c_∞ = Δc
    Δc/(K√(4Dt)) = 1/√π exp(-η^2) - η erfc(η)
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
ax[1].set_ylabel('$\mathregular{\\frac{z_+ \\nu_+ F D}{(1 - t_+^0)I} \\frac{Δc}{\sqrt{4Dt}}}$', fontsize=20)



# In[2]:
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/')

data_TP_model = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])
print(data_TP_model.keys())


# In[3]:
'''
    Comparison between
        1. analytical solution whihch assumes an infinite diffusion length, i.e. c(t,x=∞) = c_bulk
            (solid colored lines)
        2. Numerical simulation of a symmetric cell, i.e. diffusion length is not ∞

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
K = (1 - trans_Li)*i_applied_cm2/(z_Li * nu_Li * Fconst * D)


for time in [0, 0.1, 10, 50, 180]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[1].plot(data_['Position'], data_['Delta_Conc'])
    ax[1].plot(data_['Position'], func(K, D, time, data_['Position']), '--k', linewidth = 1.5)


for time in [180, 190, 230, 180*2, data_TP_model['Time'].max()]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[2].plot(data_['Position'], data_['Delta_Conc'])
    ax[2].plot(data_['Position'], func(K, D, time, data_['Position'], ti=180), '--k', linewidth = 1.5)

ax[1].set_xlim([-1e-2, 0.2])
ax[2].set_xlim([-1e-2, 0.2])



# x = 0 data
data_ = data_TP_model[data_TP_model['Position'] == 0.0]
data_1 = data_[data_['Time'] <= 180]
data_2 = data_[data_['Time'] >= 180]
ax[3].plot(data_1['Time'], data_1['Delta_Conc'])
ax[4].plot(data_2['Time'], data_2['Delta_Conc'])

ylim = ax[3].get_ylim()
ax[4].set_ylim(ylim)


ax[3].plot(data_1['Time'], func_x0(K, D, data_1['Time']), '--k', linewidth = 1.5)
ax[4].plot(data_2['Time'], func_x0(K, D, data_2['Time'], ti=180), '--k', linewidth = 1.5)

ax[1].set_ylabel('$\mathregular{\Delta c}$', fontsize=20)
ax[3].set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=20)
ax[2].set_yticklabels([])
ax[4].set_yticklabels([])

ax[1].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[2].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[3].set_xlabel('$t$', fontsize=20)
ax[4].set_xlabel('$t$', fontsize=20)

fig.tight_layout()

#



# In[3]:
rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

for time in [0, 0.1, 10, 50, 180]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[1].plot(data_['Position'], data_['Delta_Conc'])


for time in [180, 190, 230, 180*2, data_TP_model['Time'].max()]:
    data_ = data_TP_model[data_TP_model['Time'] == time]

    ax[2].plot(data_['Position'], data_['Delta_Conc'])

ax[1].set_xlim([-2e-2, 0.5])
ax[2].set_xlim([-2e-2, 0.5])

data_ = data_TP_model[data_TP_model['Position'] == 0.0]
data_1 = data_[data_['Time'] <= 180]
data_2 = data_[data_['Time'] >= 180]
ax[3].plot(data_1['Time'], data_1['Delta_Conc'])
ax[4].plot(data_2['Time'], data_2['Delta_Conc'])

ylim = ax[3].get_ylim()
ax[4].set_ylim(ylim)

ax[1].set_ylabel('$\mathregular{\Delta c}$', fontsize=20)
ax[3].set_ylabel('$\mathregular{\Delta c \ (x = 0)}$', fontsize=20)
_=ax[2].set_yticklabels([])
_=ax[4].set_yticklabels([])

ax[1].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[2].set_xlabel('$\mathregular{x}$', fontsize=20)
ax[3].set_xlabel('$t$', fontsize=20)
ax[4].set_xlabel('$t$', fontsize=20)

fig.tight_layout()

# In[5]:
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/Symmetric/')

data_TP_model_sym = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])
print(data_TP_model_sym.keys())

edge_ = data_TP_model_sym['Position'] == max(data_TP_model_sym['Position'])
data_TP_model_sym_edge = data_TP_model_sym[edge_]

new_index = range(len(data_TP_model_sym_edge))
data_TP_model_sym_edge = data_TP_model_sym_edge.reset_index(drop=True)

Pot_NJ = data_TP_model_sym_edge['Potential']
time = data_TP_model_sym_edge['Time'] / 3600.

ax, fig = axes(1, rows=2, columns=3)
# #
index = data_TP_model_sym_edge.index
relaxation = data_TP_model_sym_edge['Current'] == 0.0
relax_indices = index[relaxation]
relax_indices = relax_indices.tolist()

rest_list = []
for k, g in groupby(enumerate(relax_indices), lambda ix : ix[0] - ix[1]):
    rest_list.append(list(map(itemgetter(1), g)))

# #
ax[1].plot(time, abs(Pot_NJ)*1e3, 'b-', markersize=1)
ax[1].plot(time[relaxation], abs(Pot_NJ)[relaxation]*1e3, 'r.', markersize=1)

ax[2].plot(time[relaxation], abs(Pot_NJ)[relaxation]*1e3, 'r.', markersize=1)
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

    ax[5].plot(tau, Pot, '-')

    fit_window = (tau > 0.4) & (tau < 0.6)
    p = np.polyfit(tau[fit_window], Pot[fit_window], 1)

    tau_hat = np.linspace(0,1)
    ax[5].plot(tau_hat, np.polyval(p, tau_hat), 'k--', linewidth = 1.5)

    del_U = np.polyval(p, 1) - np.polyval(p, 0)
    ax[6].plot(I_app*np.sqrt(t_i)*1e6, del_U, 'ks')

    I_sqrt_ti.append(I_app*np.sqrt(t_i)*1e6)
    del_U_list.append(del_U)

#
# print(x.values.astype(int))
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
ax[6].set_xlabel('$\mathregular{I t_i^{1/2} \ (μA \ cm^{-2} \ s^{1/2})}$')
ax[6].set_xlim(xmin=0)
ax[6].set_ylim(ymin=0)

p = np.polyfit(I_sqrt_ti, del_U_list, 1)
I_hat = np.linspace(0, max(I_sqrt_ti))
ax[6].plot(I_hat, np.polyval(p, I_hat), 'k--', linewidth=1)

fig.tight_layout()


# In[5]:
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/Symmetric/')

data_TP_model_sym = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

for time in [0, 0.1, 10, 50, 180]:
    data_ = data_TP_model_sym[data_TP_model_sym['Time'] == time]

    ax[1].plot(data_['Position'], data_['Delta_Conc'])
    ax[1].plot(data_['Position'], func(K, D, time, data_['Position'], ti=180), '--k', linewidth = 1.5)


for time in [180, 190, 230, 180*2, data_TP_model_sym['Time'].max()]:
    data_ = data_TP_model_sym[data_TP_model_sym['Time'] == time]

    ax[2].plot(data_['Position'], data_['Delta_Conc'])
    ax[2].plot(data_['Position'], func(K, D, time, data_['Position'], ti=180), '--k', linewidth = 1.5)

# ax[1].set_xlim([-5e-2, 1+5e-2])
# ax[2].set_xlim([-5e-2, 1+5e-2])
ax[1].set_xlim([-2e-2, 0.25])
ax[2].set_xlim([-2e-2, 0.25])
ax[1].set_ylim(ylim)
ax[2].set_ylim(ylim)

data_ = data_TP_model_sym[data_TP_model_sym['Position'] == 0.0]
data_1 = data_[data_['Time'] <= 180]
data_2 = data_[data_['Time'] >= 180]
ax[3].plot(data_1['Time'], data_1['Delta_Conc'])
ax[4].plot(data_2['Time'], data_2['Delta_Conc'])

ylim = ax[3].get_ylim()
ax[4].set_ylim(ylim)

ax[3].plot(data_1['Time'], func_x0(K, D, data_1['Time'], ti=180), '--k', linewidth = 1.5)
ax[4].plot(data_2['Time'], func_x0(K, D, data_2['Time'], ti=180), '--k', linewidth = 1.5)


# In[6]:
os.chdir('/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/18_TransferencePolarization/Symmetric/')

data_TP_model_sym = pd.read_table('Time_Conc_Position.txt', delim_whitespace=True, header=[0], skiprows=[1])

rows, columns = [2, 2]
ax, fig = axes(1, rows, columns)

D = 1e-6
trans_Li = 0.2
z_Li = +1
nu_Li = 1
i_applied_cm2 = 3e-3
K = (1 - trans_Li)*i_applied_cm2/(z_Li * nu_Li * Fconst * D)

for time in [0, 0.1, 10, 50, 180]:
    data_ = data_TP_model_sym[data_TP_model_sym['Time'] == time]

    ax[1].plot(data_['Position'], data_['Delta_Conc'])
    ax[1].plot(data_['Position'], func(K, D, time, data_['Position'], ti=180), '--k', linewidth = 1.5)


for time in [180, 190, 230, 180*2, data_TP_model_sym['Time'].max()]:
    data_ = data_TP_model_sym[data_TP_model_sym['Time'] == time]

    ax[2].plot(data_['Position'], data_['Delta_Conc'])
    ax[2].plot(data_['Position'], func(K, D, time, data_['Position'], ti=180), '--k', linewidth = 1.5)

# ax[1].set_xlim([-5e-2, 1+5e-2])
# ax[2].set_xlim([-5e-2, 1+5e-2])
ax[1].set_xlim([-2e-2, 1])
ax[2].set_xlim([-5e-3, 0.1])
# ax[1].set_ylim(ylim)
ax[2].set_ylim(ylim)

data_ = data_TP_model_sym[data_TP_model_sym['Position'] == 0.0]
data_1 = data_[data_['Time'] <= 180]
data_2 = data_[data_['Time'] >= 180]
ax[3].plot(data_1['Time'], data_1['Delta_Conc'])
ax[4].plot(data_2['Time'], data_2['Delta_Conc'])

ylim = ax[3].get_ylim()
ax[4].set_ylim(ylim)

ax[3].plot(data_1['Time'], func_x0(K, D, data_1['Time'], ti=180), '--k', linewidth = 1.5)
ax[4].plot(data_2['Time'], func_x0(K, D, data_2['Time'], ti=180), '--k', linewidth = 1.5)


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
