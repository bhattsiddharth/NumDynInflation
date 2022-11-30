#########################################################################################################
#########################################################################################################
#
# This script is for plotting the phase space behaviour of the inflaton for a given model
# Please refer to <arXiv link> for explaination of variables and instructions for using the code
#
#########################################################################################################
#########################################################################################################


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from matplotlib import rcParams
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator



#########################################################################################################
# The plot settings are defined in this section
#
# It is recommended that you don't alter this section unless you're familiar with MatPlotLib and LaTeX
#########################################################################################################

# activate LaTeX text rendering
plt.rc('text', usetex=True)

#LaTex settings
plt.rcParams['text.latex.preamble']=r'\usepackage{amsmath}'
plt.rcParams['text.latex.preamble'] = r'\boldmath'


###Plot settings:

# dimensions of figure
plt.rcParams['figure.figsize'] = (10, 7)

# font
#plt.rcParams['font.size'] = 10
#plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
plt.rcParams['axes.titlesize'] = 1.4*plt.rcParams['font.size']
plt.rcParams['legend.fontsize'] = 1.4*plt.rcParams['font.size']
plt.rcParams['xtick.labelsize'] = 1.4*plt.rcParams['font.size']
plt.rcParams['ytick.labelsize'] = 1.4*plt.rcParams['font.size']

# dots per inch: dpi
#plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']

# tick sizes
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.width'] = 1

fig = plt.figure()
ax = fig.add_subplot(111)

#ticks position setting
#plt.gca().xaxis.set_ticks_position('bottom')
#plt.gca().yaxis.set_ticks_position('left')
#ax.tick_params(labeltop=False, labelright=True)

#legends
plt.rcParams['legend.frameon'] = True
#plt.rcParams['legend.loc'] = 'center left'

# width of axes
plt.rcParams['axes.linewidth'] = 1

#border setting
#plt.gca().spines['right'].set_color('none')
#plt.gca().spines['top'].set_color('none')

# gridlines
plt.grid(which='major', color='lightgrey', linestyle='-', linewidth=1, zorder=0)




#########################################################################################################
# The model of inflation is defined in this section
#########################################################################################################

# This term defines one unit of time [ T = t * m_p * S ] where t is the actual cosmic time
S = 5e-5 


# parameters used in the potential function
M = 5.9e-6 
v0 = 0.5*M**2


# dimensionless potential function and its derivatives
def f(x):
    return x**2

def dfdx(x):
    return 2*x

def d2fdx2(x):
    return 2



#########################################################################################################
# In this section we set the various initial conditions for inflation and define our dynamical equations 
# We solve the equations using the function scipy.integrate.odeint 
# After getting a numerical solution for the dynamical quantities, we plot the phase space
#########################################################################################################


### The dynamical variables are defined as follows:
#
# x : dimensionless field value [ \phi / m_p ]
# y : dimensionless field velocity [ dx/dT or \dot\phi / (m_p ^2 * S) ]
# N : number of e-folds of expansion elapsed 
# z : dimensionless hubble parameter [ H / (S * m_p) ]


# the system of differential equations to be solved
def sys(var, T):
    [x, y, z, N] = var

    # Note that all derivatives are taken wrt the scaled, dimenstionless cosmic time T

    dxdT = y
    dydT = -3*z*y - v0*dfdx(x)/S**2 
    dzdT = -0.5*y**2 #-z**2 + (v0*f(x)/S**2 - y**2)/3 # 
    dNdT = z
    
    return [dxdT, dydT, dzdT, dNdT]

# initial value of the Hubble parameter / initial energy scale 
zi = 3e-3/S


### The initial conditions are varied as follows:
# The initial Hubble parameter value is kept fixed for all cases
# We select different values of xi and find the corresponding value of yi

T = np.linspace(0, 1000, 100000)

for j in [-1,1]:
    for xi in np.arange(-10,12,2):
        yi = j*np.sqrt(6*(zi**2 - (v0*f(xi)/(3*S**2))))
        Ni = 0

        sol = odeint(sys, [xi,yi,zi,Ni], T, rtol=3e-14, atol=2e-35, mxstep=900000000)
        x, y, z, N = np.transpose(sol)
        phi, vphi, H = x, y*S, z*S

        print('%.1f\t%.4f\t%.4f'%(xi,yi,zi))
        plt.plot(phi, vphi, 'k', lw=2)

T = np.linspace(0, 20000, 100000)

for j in [-1,1]:
    yi = 0
    xi = j*np.sqrt(3 * zi**2 * S**2 / v0)
    Ni = 0

    sol = odeint(sys, [xi,yi,zi,Ni], T, rtol=3e-14, atol=2e-35, mxstep=900000000)
    x, y, z, N = np.transpose(sol)
    phi, vphi, H = x, y*S, z*S

    print('%.1f\t%.4f\t%.4f'%(xi,yi,zi))
    plt.plot(phi, vphi, 'g', lw=4)

plt.axvline(0, color='grey')
plt.axhline(0, color='grey')

### After you obtain a plot, you can manually set the limits on x and y axes to display the portion of the plot you are 

plt.xlabel(r"$\phi/m_p$", fontsize = 22) 
plt.ylabel(r"$\dot\phi/m_p^2$", fontsize = 24)

plt.show()
plt.show()