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

figure_dir = '/Users/nicholasbrady//Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/'

colors = ['blue', 'darkorange', 'green', 'purple']
markers = ['o', 's', '^']

# In[1]:
rows, columns = [1, 3]
ax, fig = axes(1, rows=rows, columns=columns)

impedance_data_file_list = [ 'Impedance_Data_Reflective.txt',
                             'Impedance_Data_Transmissive.txt',
                             'Impedance_Data_Symmetric.txt']

for a, file_ in enumerate(impedance_data_file_list):
    impedance_df = read_Fortran(file_)
    ax_ = ax[a+1]
    ax_.set_xlabel('Real(c)')
    ax_.set_ylabel('Imaginary(c)')

    for o, diff_val in enumerate([1e-2, 1e-3]):
        diff_df = impedance_df[impedance_df['Diff_Coef'] == diff_val]
        ax_.plot(diff_df['Real(c)'], abs(diff_df['Imag(c)']), 'o-', fillstyle='none', linewidth=1)

        _ = np.linspace(0, 200)
        ax_.plot(_, _, 'k--', linewidth=1)


ax[1].set_xlim(-1, 50)
ax[1].set_ylim(-1, 50)

ax[2].set_xlim(-2, 110)
ax[2].set_ylim(-2, 110)

ax[3].set_xlim(-2, 110)
ax[3].set_ylim(-2, 110)

ax[1].set_title('Reflective BC: $𝐍_{x=L} = 0$')
ax[2].set_title('Transmissive BC: $c_{x=L} = 0$')
ax[3].set_title('Symmetric BC: $𝐍_{x=L} = 𝐍_{x=0}$')

ax[1].text(2, 42, '$𝐍_{x=0} = 1$',         ha='left', va='top')
ax[1].text(2, 47, '$0 = D c\'\' - jωc$',   ha='left', va='top')

fig.savefig(os.path.join(figure_dir, 'Impedance_SimpleDiffusion.pdf'),
            format='pdf', dpi=100, bbox_inches = "tight")

# In[2]:
df_position = read_Fortran('Time_Conc_Position_Transmissive.txt')

rows, columns = [1, 2]
ax, fig = axes(1, rows=rows, columns=columns)

frequencies = np.unique(df_position['Frequency'])

# freq_df = df_position[df_position['Frequency'] == frequencies[0]]
# ax[1].plot(freq_df['Position'].values, freq_df['Real(c)'].values, 'o-', linewidth=1, fillstyle='none')

for i in range(0, 20, 4):
    freq_df = df_position[df_position['Frequency'] == frequencies[i]]
    real_c = freq_df['Real(c)'].values
    ax[1].plot(freq_df['Position'].values, real_c / max(abs(real_c)), 'o-', linewidth=1, fillstyle='none')

    real_c = abs(freq_df['Real(c)'].values)
    ax[2].semilogy(freq_df['Position'].values, real_c / max(real_c), 'o-', linewidth=1, fillstyle='none')


# print(type(freq_df['Real(c)'].values[0]))
