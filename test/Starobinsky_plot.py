#Axis labels using LaTex format
#Make labels bold with appropriate spacing
#Use of log/semilog scale

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from matplotlib import rcParams
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator

# activate latex text rendering
plt.rc('text', usetex=True)

#LaTex setting
plt.rcParams['text.latex.preamble']=r'\usepackage{amsmath}'
plt.rcParams['text.latex.preamble'] = r'\boldmath'

#Plot setting:
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plt.rcParams['figure.figsize'] = (10, 7)
#plt.rcParams['font.size'] = 10
#plt.rcParams['font.family'] = 'Times New Roman'

plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
plt.rcParams['axes.titlesize'] = 1.4*plt.rcParams['font.size']
#plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
plt.rcParams['xtick.labelsize'] = 1.4*plt.rcParams['font.size']
plt.rcParams['ytick.labelsize'] = 1.4*plt.rcParams['font.size']
# dots per inch: dpi
#plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']

plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.width'] = 1

#legends
plt.rcParams['legend.frameon'] = True
#plt.rcParams['legend.loc'] = 'center left'
plt.rcParams['legend.fontsize'] = 14*plt.rcParams['font.size']

plt.rcParams['axes.linewidth'] = 1

#border setting
#plt.gca().spines['right'].set_color('none')
#plt.gca().spines['top'].set_color('none')

#ticks position setting
#plt.gca().xaxis.set_ticks_position('bottom')
#plt.gca().yaxis.set_ticks_position('left')
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.tick_params(labeltop=False, labelright=True)

plt.grid(which='major', color='lightgrey', linestyle='-', linewidth=1, zorder=0)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



S = 5e-5 # for scaling time


# parameters used in the potential function
b = np.power(2.0/3, 0.5) # \sqrt(2/3)
v0 = 0.96*np.power(10.0, -10)

def f(x):
    return np.power(1 - np.exp(-b*x),2)

def dfdx(x):
    return 2*b*np.exp(-b*x)*(np.exp(b*x) - 1)

def d2fdx2(x):
    return -(4.0/3)*np.exp(-b*x)*(np.exp(b*x) - 2)



#   Select the sequence to execute:

switch = 0

#   background          0
#   calculate scales    1
#   scalar modes        2
#   tensor modes        3
#   power spectrum      4
#   ratio               5




if switch == 0:
#''' # background
    data = np.loadtxt('data/Starobinsky_MS_background.txt') # T,N,Ne,x,y,etaH,epsH,Ps,Pt 9
    data1 = np.loadtxt('data/Starobinsky_MS_scales.txt') # Ne,x,y,z,A,aH 6

    x = data[:,3]
    y = data[:,4]
    z = np.sqrt(y**2/6 + (v0*f(x)/(3*S**2)))
    meff = 2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2
    aH = data1[:,5]

    def i(Ne):
        return np.max(np.where(data[:,2]>=Ne))

    ''' # to calculate n_s and r
    ns = 1 + 2*data[:,5] - 4*data[:,6]
    r = 16*data[:,6]
    print('%.4f \t %.4f'%(1 + 2*data[i(50),5] - 4*data[i(50),6], 1 + 2*data[i(60),5] - 4*data[i(60),6]))
    print('%.4f \t %.4f'%(16*data[i(50),6], 16*data[i(60),6]))
    '''

    ''' # for T vs N
    plot1 = plt.plot(data[:,0],data[:,1])
    plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    plt.text(100,20, r'$\textbf{accelerated expansion}$', fontsize=20, color='purple', rotation=35, bbox=dict(fc='w', ec='k', boxstyle='round'))
    plt.text(690,35, r'$\textbf{end of inflation}$', fontsize=20, color='purple', rotation=90, bbox=dict(fc='w', ec='k', boxstyle='round'))
    plt.axhline(y=79.4135, ls=':', lw=2, color='k')
    plt.axvline(x=730, ls='--', lw=2, color='k')
    plt.text(350,25, r'$a\sim \mathrm{e}^{Ht}$', fontsize=20, bbox=dict(fc='w', ec='k', pad=10))
    plt.text(770,30, r'$\textbf{decelerated expansion}$', fontsize=20, color='purple', rotation=90, bbox=dict(fc='w', ec='k', boxstyle='round'))
    plt.xlim(0,800)
    plt.ylim(0,85)
    plt.xlabel(r'$T$', fontsize = 22) 
    plt.ylabel(r'$N$', fontsize = 24)
    '''

    ''' # for phase plot
    plot1 = plt.plot(data[:,3],data[:,4]*S)
    plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    plt.axvline(x=0.615, ls='--', lw=2, color='k')
    plt.text(0.7, -3e-6, r'$\textbf{end of inflation}$', fontsize=20, color='purple', rotation=90)
    plt.xlabel(r'$\phi/m_p$', fontsize = 22) 
    plt.ylabel(r'$\dot\phi/m_p^2$', fontsize = 24)
    '''

    ''' # for slow-roll parameters
    plot1 = plt.plot(data[:,2],np.abs(data[:,5]))
    plt.setp(plot1, color='r', linewidth=3.5, linestyle='-')

    plot1 = plt.plot(data[:,2],data[:,6])
    plt.setp(plot1, color='g', linewidth=3.5, linestyle='-')

    ax.add_patch(pch.Rectangle((0,1), 70,10, fc='brown', alpha=0.5))
    plt.text(5,2, r'$\epsilon_H,\,\eta_H > 1$', fontsize=20, bbox=dict(fc='w', ec='k'))
    plt.axhline(y=1, ls='--', color='k', lw=2)
    plt.text(0.5,2e-3, r'$\textbf{end of inflation}$', fontsize=20, color='purple', verticalalignment='center', rotation=90)
    
    plt.text(30,6e-2, r'$|\eta_H|$', color='k', fontsize=18, bbox=dict(fc='w', ec='r', boxstyle='round'))
    ax.annotate('', xy=(25, 0.0365), xytext=(23, 0.0394), arrowprops=dict(arrowstyle='<-', lw=5, color='r'))
    
    plt.text(20,3e-3, r'$\epsilon_H$', color='k', fontsize=18, bbox=dict(fc='w', ec='g', boxstyle='round'))
    ax.annotate('', xy=(30, 7.25e-4), xytext=(28, 8.26e-4), arrowprops=dict(arrowstyle='<-', lw=5, color='g'))

    plt.text(40,5e-3, r'$\textbf{slow-roll: }$', fontsize=20, color='purple')
    plt.text(37,2e-3, r'$\epsilon_H,\,|\eta_H|\ll 1$', fontsize=20, bbox=dict(fc='w', ec='k'))
    
    plt.axvline(x=60, ls='--', lw=2, color='b')
    ax.add_patch(pch.Rectangle((55,1e-5), 10,10, fc='grey', alpha=0.5))
    plt.plot([60,60],[np.abs(data[i(60),5]),data[i(60),6]], 'b*', ms=14)
    plt.text(61,2e-3, r'$\textbf{pivot scale}$', fontsize=18, color='purple', rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))
    plt.text(60,1.3e-5, r'$\textbf{CMB window}$', fontsize=12, color='purple', horizontalalignment='center', bbox=dict(fc='w', ec='k'))

    plt.yscale('log')

    plt.xlim(0.0,70.0)
    plt.ylim(1e-5,10)

    plt.xlabel(r'$N_e$', fontsize=22)
    '''

    ''' # for z''/z
    plot1 = plt.plot(data1[:,0], z**(-2) *meff, label=r'$(aH)^{^{-2}}\frac{z^{\prime\prime}}{z}$')
    plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    plot1 = plt.plot(data[:,2], data[:,7], label=r'$\eta_H$')
    plt.setp(plot1, color='b', linewidth=4, linestyle='--')

    plot1 = plt.plot(data[:,2], 2 - data[:,6], label=r'$(aH)^{^{-2}}\frac{a^{\prime\prime}}{a}$')
    plt.setp(plot1, color='r', linewidth=4, linestyle='--')

    plt.xlim(0,70)
    plt.ylim(0,3)

    plt.xlabel(r'$N_e$', fontsize = 22) 
    plt.ylabel(r'$(aH)^{-2} z^{\prime\prime}/z$', fontsize = 24)
    plt.legend(frameon=True, fontsize = 16, borderpad = 1)
    '''

    ''' # for x,y,z vs Ne
    plot1 = plt.plot(data[:,2], S*np.sqrt(data[:,4]**2/6 + (v0*f(data[:,3])/(3*S**2))), label=r'$H/m_p$')
    plt.setp(plot1, color='g', linewidth=4, linestyle='-')


    #plot1 = plt.plot(data[:,2], data[:,3], label=r'$\phi/m_p$')
    #plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    #plot1 = plt.plot(data[:,2], S*data[:,4], label=r'$\dot\phi/m_p^2$')
    #plt.setp(plot1, color='r', linewidth=4, linestyle='-')


    plt.axvline(x=60, ls='--', lw=2, color='b')
    ax.add_patch(pch.Rectangle((55,-5), 10,20, fc='grey', alpha=0.5))


    #plt.plot([60,60], [data[i(60),3],S*data[i(60),3]], 'b*', ms=14)
    #plt.text(56,3, r'$\textbf{pivot scale}$', fontsize=24, color='purple', rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))
    #plt.text(60,1, r'$\textbf{CMB window}$', fontsize=20, color='purple', horizontalalignment='center', bbox=dict(fc='w', ec='k'))


    plt.plot(60, S*np.sqrt(data[i(60),4]**2/6 + (v0*f(data[i(60),3])/(3*S**2))), 'b*', ms=14)
    plt.text(56,3e-6, r'$\textbf{pivot scale}$', fontsize=24, color='purple', rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))
    plt.text(60,0.5e-6, r'$\textbf{CMB window}$', fontsize=20, color='purple', horizontalalignment='center', bbox=dict(fc='w', ec='k'))

    ax.add_patch(pch.Rectangle((3,2e-6), 47,1.8e-6, fill=True, alpha=1, fc='w', ec='k', lw=2, zorder=10))
    plt.text(4,2e-6, r'$$H^2 = \frac{1}{3m_p^2} \left[\frac{\dot\phi^2}{2}+V(\phi)\right]$$', fontsize=32, zorder=20)#, bbox=dict(fc='w', ec='k'))

    plt.xlim(0,75)
    #plt.ylim(-1,6)
    plt.ylim(0,7e-6)
    plt.xlabel(r'$N_e$', fontsize=32)
    plt.ylabel(r'$H/m_p$', fontsize=34)
    #plt.legend(frameon=True, fontsize = 24, borderpad = 1, loc='best') #legend format
    ''' 

    ''' # for n_S and n_T
    epsH = data[:,6]
    etaH = data[:,5]

    ns = 2*etaH - 4*epsH
    nt = 2*epsH

    plot1 = plt.plot(data[:,2], ns)
    plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    #plot1 = plt.plot(data[:,2], nt)
    #plt.setp(plot1, color='g', linewidth=4, linestyle='-')

    plt.axvline(60, color='b', linewidth=2, linestyle='--')
    plt.plot(60.01,0.9674-1, 'b*', ms=17) # ns
    #plt.plot(60.01,0.0004, 'b*', ms=16) # nt
    ax.add_patch(pch.Rectangle((55,-0.3), 10,0.7, fc='grey', alpha=0.5))

    #plt.text(55,4.11, r'$\times 10^{-3}$', fontsize=22, horizontalalignment='left', verticalalignment='top')

    #plt.yscale('log')
    #plt.grid(which='both', color='lightgrey', linestyle='-', linewidth=1)
    plt.xlim(0,70)
    plt.ylim(-0.3,0)
    #plt.ylim(1e-5,1)
    plt.xlabel(r'$N_e$', fontsize=34)
    plt.ylabel(r'$n_{_S}-1$', fontsize=34)
    #plt.ylabel(r'$|n_{_T}|$', fontsize=36)
    '''
#'''



if switch == 1:
#''' ### to calculate scales
    data = np.loadtxt('data/Starobinsky_MS_scales.txt') # Ne,x,y,z,A,aH 6

    def i(Ne):
        return np.max(np.where(data[:,0]>=Ne))

    #for N in range(0,76):
    #    print('%d\t%3g\t%3g\t%3g'%(data[i(N),0],data[i(N),5],data[i(N),1],data[i(N),2]))

    plot1=plt.plot(data[:,0], data[:,5]) 
    plt.setp(plot1, color='r', linewidth=2.5, linestyle='-')

    #plt.xscale('log')
    plt.yscale('log')

    #plt.xlim(0,70)
    #plt.ylim(2.0*10*(-10),2.3*10*(-9),)

    plt.xlabel(r'$N_e$', fontsize = 22) 
    plt.ylabel(r'$aH$', fontsize = 24)
#'''



if switch == 2:
#''' ### to calculate P_S at k
    for Nk in [60]: #np.arange(0,72,2):
        file_name = 'data/modes/Starobinsky_MS_amplitude%.1f.txt'%Nk
        data = np.loadtxt(file_name) #Ne,aHk,P_S,P_T 4

        def i(Ne):
            return np.max(np.where(data[:,0]>=Ne))

        def j(aHk):
            return np.max(np.where(data[:,1]<=aHk))

        print('%.1f\t%3g\t%3g'%(Nk,data[i(Nk-5),2],data[i(Nk-5),3]))

        c='k'
        if Nk in np.arange(48,72,2): c='b'
        if Nk in np.arange(24,48,2): c='g'
        if Nk in np.arange(0,24,2): c='r'

        plt.axhline(y=2.105e-9, ls=':', lw=2, color='b')

        plot1 = plt.plot(data[:,1], data[:,2], label='Ne = %.1f'%Nk) # x axis: aH/k
        plt.setp(plot1, color='g', linewidth=4, linestyle='-')

        #plot1 = plt.plot(data[:,0] - 79.4135 + Nk+5, data[:,2], label='Ne = %.1f'%Nk)
        #plt.setp(plot1, color=c, linewidth=2.5, linestyle='-')

        #plt.plot(data[j(1.0),0]- 79.4135 + Nk+5, data[j(1.0),2], marker='o', markersize=5, color='k')

        ax.annotate('', xy=(0.12, 1.55e-7), xytext=(0.103, 2.1e-7), arrowprops=dict(arrowstyle='->', lw=5, color='g'))

        plt.axvline(x=1, ls='--', lw=2.5, color='k')

        plt.text(1.5, 1e-6, r'$\textbf{Hubble exit}$', fontsize=18)
        plt.text(1.6, 5e-7, r'$k=aH$', fontsize=18, bbox=dict(facecolor='none', edgecolor='black'))

        plt.text(9e-3, 5e-9, r'$\textbf{sub-Hubble}$', fontsize=20, color='purple')
        plt.text(20, 5e-9, r'$\textbf{super-Hubble}$', fontsize=20, color='purple')

        #plt.text(0.1, 2e-5, r'$\textbf{mode corresponds to CMB scale }(k_*)$', fontsize=20, color='b', bbox=dict(facecolor='w', edgecolor='b', boxstyle='round'))

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim(7e-3,5e2)
    #plt.ylim(1e-12,1e-8)

    plt.xlabel(r'$aH/k$', fontsize = 22) 
    plt.ylabel(r'$k^3\left|\zeta_k\right|^2/2\pi^2$', fontsize = 24) 
#'''



if switch == 3:
#''' ### to calculate P_T at k
    for Nk in [60]: #np.arange(70,-1,-2):
        file_name = 'data/modes/Starobinsky_MS_amplitude%.1f.txt'%Nk
        data = np.loadtxt(file_name) #Ne,aHk,P_S,P_T 4

        def i(Ne):
            return np.max(np.where(data[:,0]>=Ne))

        def j(aHk):
            return np.max(np.where(data[:,1]<=aHk))

        #print('%.1f\t%3g\t%3g'%(Nk,data[i(Nk-2),2],data[i(Nk-2),3]))

        c='k'
        if Nk in np.arange(48,72,2): c='b'
        if Nk in np.arange(24,48,2): c='g'
        if Nk in np.arange(0,24,2): c='r'

        plot1 = plt.plot(data[:,0] - 79.41356 + Nk+5, data[:,3], label='Ne = %.1f'%Nk)
        plt.setp(plot1, color=c, linewidth=2, linestyle='-')

        plt.plot(data[j(1.0),0]- 79.41356 + Nk+5, data[j(1.0),3], marker='o', markersize=5, color='k')

    #plt.text(25, 2e-11, r'$\textbf{Power spectrum spans almost }$', fontsize=18)
    #plt.text(25, 9e-12, r'$\textbf{two orders of magnitude }$', fontsize=18)
    
        
    #plt.xscale('log')
    plt.yscale('log')

    #plt.xlim(0,75)
    #plt.ylim(2e-12,2e-9)

    plt.xlabel(r'$N_e$', fontsize = 22) 
    plt.ylabel(r'$4k^3\left|h_k\right|^2/\pi^2$', fontsize = 24) 
#'''



if switch == 4:
#''' ### To plot Power spectrum vs Ne
    data = np.loadtxt('data/Starobinsky_MS_power.txt') # Ne,P_S,P_T 3
    data1 = np.loadtxt('data/Starobinsky_MS_background.txt') # T,N,Ne,x,y,etaH,epsH,Ps,Pt 9

    def i(Ne):
        return np.max(np.where(data1[:,0]<=Ne))

    plt.plot(data[:,0],data[:,1], 'go', ms=7, label='Scalar (numerical)')
    plt.plot(data1[:,2],data1[:,7], 'g-',lw=3, label='Scalar (slow-roll)')

    plt.plot(data[:,0],data[:,2], 'ro', ms=7, label='Tensor (numerical)')
    plt.plot(data1[:,2],data1[:,8], 'r-',lw=3, label='Tensor (slow-roll)')

    plt.axhline(y=2.1e-9, ls=':', lw=2, color='b')
    plt.plot([60,60],[2.1e-9, data1[i(60),8]], 'b*', ms=14)
    plt.text(30,2.5e-9, r'$\textbf{CMB normalization}$', fontsize=18, color='purple')
    #plt.text(6,3e-9, r'$A_{_S}=2.1\times 10^{-9}$', fontsize=14, color='k', bbox=dict(fc='w', ec='k'))

    plt.text(58,1.2e-9, r'$n_{_S} \simeq 0.967$', fontsize=14, bbox=dict(fc='w', ec='k'), horizontalalignment='left')
    plt.text(58,7e-10, r'$r \simeq 0.003$', fontsize=14, bbox=dict(fc='w', ec='k'), horizontalalignment='left')

    #plt.text(10,2e-11, r'$\mathcal{P}_{_S}=\frac{1}{8\pi^2}\left(\frac{H}{m_p}\right)^2\frac{1}{\epsilon_H}$', color='k', fontsize=22, bbox=dict(lw=2 ,fc='w', ec='g', boxstyle='round'))
    #plt.text(10,2e-12, r'$\mathcal{P}_{_T}=\frac{2}{\pi^2}\left(\frac{H}{m_p}\right)^2$', color='k', fontsize=22, bbox=dict(lw=2 ,fc='w', ec='r', boxstyle='round'))

    plt.axvline(x=60, ls='--', lw=2, color='b')
    ax.add_patch(pch.Rectangle((55,1e-12), 10,(1e-8-1e-12), fc='grey', alpha=0.5))
    plt.text(57,5e-11, r'$k_*=0.05\,\mathrm{Mpc}^{-1}$', fontsize=14, rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))    
    plt.text(61,5e-11, r'$\textbf{pivot scale}$', fontsize=18, color='purple', rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))
    plt.text(60,1.3e-12, r'$\textbf{CMB window}$', fontsize=12, color='purple', horizontalalignment='center', bbox=dict(fc='w', ec='k'))

    #ax.add_patch(pch.Rectangle((65,1e-12), 10,(1e-8-1e-12), fc='brown', ec='brown', alpha=0.5))
    #plt.text(70,1e-10, r'$\textbf{large-scale fluctuations outside the horizon}$', fontsize=18, color='purple', rotation=90, verticalalignment='center', bbox=dict(fc='w', ec='k'))
    
    #ax.annotate('', xy=(38,1e-12), xytext=(35,1e-12), arrowprops=dict(arrowstyle='<-', lw=1, color='k'))

    plt.yscale('log')

    plt.xlim(0,75)
    plt.ylim(1e-12,1e-8)

    def k(Ne):
        return 0.05*np.exp(60 - Ne)

    def Ne(k):
        return 60 - np.log(k/0.05)

    #secax = ax.secondary_xaxis('top', functions=(k, Ne))
    #secax.set_xscale('log')
    #secax.set_xlabel(r'$k/\mathrm{Mpc}^{-1}$', fontsize = 22)

    plt.xlabel(r'$N_e$', fontsize = 22) 
    plt.ylabel(r'$\mathcal{P}_{S,T} = \frac{k^3}{2\pi^2} \left| \zeta_k \right|^2 _{S,T}$', fontsize = 24) 
    plt.title(r'$\textbf{Power spectra for Starobinsky Model}$', fontsize = 22)
    plt.legend(frameon=True, fontsize = 14, borderpad = 1.2) #legend format
#'''


if switch == 5:
#''' ### ratio
    data = np.loadtxt('data/Starobinsky_MS_power.txt') # Ne,P_S,P_T
    data1 = np.loadtxt('data/Starobinsky_MS_background.txt') # T,N,Ne,x,y,etaH,epsH,As,Ps,At,Pt 11

    def i(Ne):
        return np.max(np.where(data1[:,2]>=Ne))

    ratio = data[:,2]/data[:,1]

    r = CubicSpline(data[:,0], ratio)

    plot1 = plt.plot(data[:,0], ratio, label=r'$r$')
    plt.setp(plot1, color='g', linewidth=2.5)

    plot1 = plt.plot(data1[:,2], 16*data1[:,6], label=r'$16\epsilon_H$')
    plt.setp(plot1, color='r', linewidth=2.5)

    plt.plot(60, r(60), marker='o', ms=9, color='k')
    plt.text(50,1, 'r(Ne=60) = %.3g'%r(60), fontsize=16)

    plt.yscale('log')

    plt.xlim(0,70)
    plt.ylim(1e-3,1e1)

    plt.xlabel(r'$N_e$', fontsize = 22) 
    plt.ylabel(r'$\mathrm{tensor-to-scalar}$', fontsize = 24) 
#'''

plt.show()
