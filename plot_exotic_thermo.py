
# coding: utf-8

# In[ ]:

# import necessary modules
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from matplotlib import rc
import matplotlib.patches as patches

from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator
from math import floor
from mpl_toolkits.axes_grid1 import make_axes_locatable

import time
start_time = time.time()

# In[ ]:

# esthetic definitions for the plots

#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)
# matplotlib.rc('font', **font)
#matplotlib.mathtext.rcParams['legend.fontsize']='medium'
#plt.rcParams["figure.figsize"] = [8.0,6.0]
# norm = matplotlib.colors.Normalize(vmin=min(weights), vmax=max(weights))

cmap = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.Blues)

###planck 2018: problem, only the large l are binned. lowell unbinned?
lTT,DlTT_mean,DlTT_error_minus,DlTT_error_plus,DlTT_bestfit= np.loadtxt("error_Planck/Planck2018_errorTT.txt",unpack=True)
lEE,DlEE_mean,DlEE_error_minus,DlEE_error_plus,DlEE_bestfit= np.loadtxt("error_Planck/Planck2018_errorEE.txt",unpack=True)
lTE,DlTE_mean,DlTE_error_minus,DlTE_error_plus,DlTE_bestfit= np.loadtxt("error_Planck/Planck2018_errorTE.txt",unpack=True)

##2015 error
# lTT,lminTT,lmaxTT,DlTT,DlTT_error= np.loadtxt("error_Planck/Planck2015_errorTT.txt",unpack=True)
# lEE,lminEE,lmaxEE,DlEE,DlEE_error= np.loadtxt("error_Planck/Planck2015_errorEE.txt",unpack=True)
# lTTlowl,DlTTlowl,DlTT_errorup_lowl,DlTT_errordown_lowl= np.loadtxt("error_Planck/Planck2015_errorTT_lowl.txt",unpack=True)
# lEElowl,DlEElowl,DlEE_errorup_lowl,DlEE_errordown_lowl= np.loadtxt("error_Planck/Planck2015_errorEE_lowl.txt",unpack=True)

# lTT_all = np.append(lTTlowl,lTT)
# lEE_all = np.append(lEElowl,lEE)
# fract_TT_down = np.append((DlTT_errordown_lowl)/DlTTlowl,(DlTT_error)/DlTT)
# fract_TT_up = np.append((DlTT_errorup_lowl)/DlTTlowl,(DlTT_error)/DlTT)
# fract_EE_down = np.append((DlEE_errordown_lowl)/DlEElowl,(DlEE_error)/DlEE)
# fract_EE_up = np.append((DlEE_errorup_lowl)/DlEElowl,(DlEE_error)/DlEE)
# print(fract_TT_down,fract_TT_up)
# print(lTT_all,lEE_all)



#ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=DlEE_error_plus/DlEE_mean, fmt='.',color='r')



##create plot
ax_1 = plt.subplot(411)
ax_2 = plt.subplot(412, sharex = ax_1)
ax_3 = plt.subplot(413)
ax_4 = plt.subplot(414)
#plt.subplots_adjust(hspace=0)
# ax_1.set_ylim([-0.12,0.12])
# ax_2.set_ylim([-0.12,0.12])
# ax_1.set_ylim([-0.06,0.06])
# ax_2.set_ylim([-0.06,0.06])
ax_1.set_xlim([2,2500])
ax_2.set_xlim([2,2500])
# ax_1.set_ylim([30,30000])
# ax_2.set_ylim([-1,9])
# ax_2.set_ylim([-2,2])

# PUT PLANCK 2015 TT, TE, EE+lowP+lensing best-fit values
#for omega_b, n_s, A_s, tau_reio, Omega_ini_dcdm2 and H0

##LCDM BESTFIT Planck 2018####
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   'l_max_scalars':2600,
                   'n_s':0.9652,
                   'ln10^{10}A_s':3.043,
                   # 'tau_reio':0.0540,
                   'omega_b':0.02233,
                   'h':0.6737,
                   # '100*theta_s':1.041,
                   'lensing':'yes',
                   'recfast_Nz0':100000,
                   'recombination':'hyrec'
                   # 'reio_parametrization':'reio_none'
                   # 'temperature contributions':'lisw,eisw'
                   }

M = Class()
M.set(common_settings)
M.compute()
clM = M.lensed_cl(2600)
thermo_ref = M.get_thermodynamics()
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
kvec = np.logspace(-4,0,100)
T_cmb = 2.7225 #we change units for Planck
fTT_ref = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
fEE_ref = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
M.struct_cleanup()
M.empty()
timeafterref=time.time()
# z_table_to_change = [700,800,900,1000,1100,1200,1300,1400,1500]
# fractional_change_xe = 1,1,1,1,1,1,1,1,1
size_z_table_to_change = 10
zmin = 600
zmax = 1200
ax_3.set_xlim([zmin-10,zmax+10])
ax_4.set_xlim([zmin-10,zmax+10])

z_table_to_change = np.linspace(zmin,zmax,size_z_table_to_change)
# z_table_to_change = [800,900,750,900,850,900,850,1250,800,1300]
zrange = np.arange(0,2500)
xe_ref = interp1d(thermo_ref['z'],thermo_ref['x_e'])
# kappaprime_ref = interp1d(thermo_ref['z'],thermo_ref['kappa\' [Mpc^-1]'])
kappaprime_ref = interp1d(thermo_ref['z'],thermo_ref['g [Mpc^-1]'])
print(thermo_ref['z'])
# print(z_table_to_change)
plus_ou_moins = [-0.1,+0.1,1.03,0.97,1.05,0.95]
# ax_3.set_ylim([(plus_ou_moins[1]-1)*2,(plus_ou_moins[0]-1)*2])
color=['red','red','blue','blue','green','green']
###choose the value of Gamma and the number of bins in perts
for ii in np.arange(0,2):
    for i in np.arange(0,size_z_table_to_change-1):
    # for i in [0,2,4]:
    # for m_ncdm in [0.9/3]:
            print('%.2f'%(plus_ou_moins[ii]))
            if ii == 0:
                if i == 0:
                    var_color = 'r'
                else:
                    var_color = plt.cm.Reds(0.9*i/(size_z_table_to_change)+0.3)
            if ii ==1:
                if i == 0:
                    var_color = 'g'
                else:
                    var_color = plt.cm.Greens(0.9*i/(size_z_table_to_change)+0.3)
            print("~~~~~i %d computing z = [%f,%f]~~~~~"%(i,z_table_to_change[i],z_table_to_change[i+1]))
            M.set(common_settings)

            # print("here",j,plus_ou_moins[j])
            print('Here')
            M.set({
            'z_table_to_change': '%f'%z_table_to_change[i],
            'sigma_z_table_to_change':'10',
            # 'z_table_to_change': '%f,%f'%(z_table_to_change[i],z_table_to_change[i+1]),
            'fractional_change_xe': '%.2f'%(plus_ou_moins[ii]),
            # 'z_table_to_change': '800,840,843,845,847,850',
            # 'fractional_change_xe': '1.01,1.006,1.006,1.003,1.003,1',
            'size_z_table_to_change':1
            })
            M.compute()
            clM = M.lensed_cl(2600)
            ll_LCDM = clM['ell'][2:]
            clTT_LCDM = clM['tt'][2:]
            clEE_LCDM = clM['ee'][2:]
            kvec = np.logspace(-4,0,100)
            T_cmb = 2.7225 #we change units for Planck
            fTT_ourcode = interp1d(ll_LCDM,clTT_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
            fEE_ourcode = interp1d(ll_LCDM,clEE_LCDM*(ll_LCDM)*(ll_LCDM+1)/2/np.pi*(T_cmb*1.e6)**2)
            thermo_ourcode = M.get_thermodynamics()

            xe_ourcode = interp1d(thermo_ourcode['z'],thermo_ourcode['x_e'])
            # kappaprime_ourcode = interp1d(thermo_ourcode['z'],thermo_ourcode['kappa\' [Mpc^-1]'])
            kappaprime_ourcode = interp1d(thermo_ourcode['z'],thermo_ourcode['g [Mpc^-1]'])
            # print(thermo_ourcode['z'],thermo_ourcode['x_e'])
            M.struct_cleanup()
            M.empty()




            ax_1.semilogx(ll_LCDM,fTT_ourcode(ll_LCDM)/fTT_ref(ll_LCDM)-1,color=var_color)
            ax_2.semilogx(ll_LCDM,fEE_ourcode(ll_LCDM)/fEE_ref(ll_LCDM)-1,color=var_color)
            ax_3.plot(zrange,(xe_ourcode(zrange)/xe_ref(zrange)-1),color=var_color, label=r'$z=[%.0f,%.0f]$'%(z_table_to_change[i],z_table_to_change[i+1]))
            ax_4.plot(zrange,(kappaprime_ourcode(zrange)/kappaprime_ref(zrange)-1),color=var_color)







print("~~~~~compute binned cosmic variance based on ref~~~~~")

def binned_cosmic_variance (result,l_ini,width):
    central_l = l_ini+width/2
    weight_total = 0
    result = 0
    Clb = 0
    for i in range(0,int(width)):
        weight_total += (l_ini+i)*(l_ini+i+1)
        result += 2/(2*(l_ini+float(i))+1)*(l_ini+float(i))*(l_ini+float(i)+1)*(l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)*fTT_ref(l_ini+i)
        Clb += (l_ini+float(i))*(l_ini+float(i)+1)*fTT_ref(l_ini+i)
    return np.sqrt(result)/Clb

#
l_min = 2.;
l_max = 2000;
n_step = 100.;
j=0.
step = l_min
width= 30
while step < l_max:
        result = 0.0
        if step < 29:
            width = 1
            step = l_min+j*width
            j+=1
            if step == 29:
                j = 0
                l_min = 30
        else:
            width = 30
            step = l_min+j*width
            j+=1

        ax_1.add_patch(patches.Rectangle((int(step), -1*binned_cosmic_variance(result,int(step),width)), width, 2*binned_cosmic_variance(result,int(step),width),color='r', alpha=0.05))
        ax_2.add_patch(patches.Rectangle((int(step), -1*binned_cosmic_variance(result,int(step),width)), width, 2*binned_cosmic_variance(result,int(step),width),color='r',alpha=0.05))

print("~~~~~ready to plot~~~~~")

print("~~~~~print planck error bar around mean, i.e, residuals are 0 for simplicity~~~~~")


ax_1.errorbar(lTT, DlTT_mean/DlTT_mean-1, yerr=[DlTT_error_minus/DlTT_mean,DlTT_error_plus/DlTT_mean], fmt='.',color='lightgray')
ax_2.errorbar(lEE, DlEE_mean/DlEE_mean-1, yerr=[DlEE_error_minus/DlEE_mean,DlEE_error_plus/DlEE_mean], fmt='.',color='lightgray')

l_cosmic_variance = np.linspace(0,48,1000)

l_cosmic_variance_1 = np.linspace(0,30,1000)

l_cosmic_variance_2 = np.linspace(30,48,2)
slope =np.array([0.13,0.0343])

# ax_1.fill_between(l_cosmic_variance_1, -0.13,0.13, color='lightgray' )
# ax_1.fill_between(l_cosmic_variance_2, -slope,slope, color='lightgray' )
# ax_1.fill_between(lTT, -(DlTT_error)/DlTT, +(DlTT_error)/DlTT, color='lightgray')
# ax_1.fill_between(lTTlowl, (DlTT_errordown_lowl-DlTTlowl)/DlTTlowl, +(DlTT_errorup_lowl-DlTTlowl)/DlTTlowl, color='lightgray')
# ax_1.fill_between(lTT_all, fract_TT_down, fract_TT_up, color='lightgray')
# ax_2.fill_between(lEE_all, fract_EE_down, fract_EE_up, color='lightgray')
# ax_1.errorbar(lTT_all,0*lTT_all,yerr=[fract_TT_down, fract_TT_up], color='lightgray',fmt='.')
# ax_2.errorbar(lEE_all,0*lEE_all,yerr=[fract_EE_down, fract_EE_up], color='lightgray',fmt='.')






ax_1.set_xlabel(r'$\mathrm{multipole} \, \, \ell$',fontsize=15)
# ax_3.set_xlabel(r'$k$',fontsize=15)
# ax_1.set_ylabel(r'$\frac{C_\ell^\mathrm{TT}(\mathrm{old})}{C_\ell^\mathrm{TT}(\mathrm{new} )} -1$',fontsize=15)
# ax_2.set_ylabel(r'$\frac{C_\ell^\mathrm{EE}(\mathrm{old})}{C_\ell^\mathrm{EE}(\mathrm{new} )} -1$',fontsize=15)
# ax_3.set_ylabel(r'$\Delta x_e/x_e$',fontsize=15)
ax_4.set_ylabel(r'$\Delta g/g$',fontsize=15)
ax_3.set_ylabel(r'$\Delta x_e/x_e$',fontsize=15)
ax_1.set_ylabel(r'$\Delta C_l(TT)/C_l(TT)$',fontsize=15)
ax_2.set_ylabel(r'$\Delta C_l(EE)/C_l(EE)$',fontsize=15)

ax_1.set_title(r'$h=%.2f, z\in[%.f,%.f]$'%(0.6737,zmin,zmax), fontdict={'fontsize':15, 'fontweight': 'medium'})# ax_3.tick_params(axis="x", labelsize=15)
ax_2.tick_params(axis="x", labelsize=15)
ax_1.tick_params(axis="x", labelsize=15)
ax_2.tick_params(axis="y", labelsize=15)
ax_1.tick_params(axis="y", labelsize=15)
ax_3.tick_params(axis="y", labelsize=15)
ax_4.tick_params(axis="y", labelsize=15)
ax_3.tick_params(axis="x", labelsize=15)
ax_4.tick_params(axis="x", labelsize=15)

#ax_3.legend(frameon=False,fontsize = 12,loc='upper left',borderaxespad=0.)

plt.show()
plt.clf()



#

# In[ ]:
