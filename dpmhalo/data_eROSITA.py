# Function adapted from Yi Zhang's data that she sent me.  

import numpy as np
import os, sys, glob
import astropy.units as u
import matplotlib.pyplot as p
from astropy.cosmology import FlatLambdaCDM

ModelDir = '.'

def Zhang_2024_eRASS4_Xray(massbinname):

    data_file = ModelDir + '/Data/profiles_feb21_forErwin/profile_CENhalo_CGM_XRB_%s.txt'%massbinname
    if(massbinname=="M31"):
        M_200m = 1.4e+12
    if(massbinname=="2M31"):
        M_200m = 4.6e+12
    if(massbinname=="4M31"):
        M_200m = 1.35e+13

    lM_200m = np.log10(M_200m)
    
    bins_mean,bins_low,bins_up,profile_gal_cen,profile_gal_cen_err,psf_gal_cen=np.transpose(np.loadtxt(data_file,unpack=True))[:6]
    profile_subsatboost=np.transpose(np.loadtxt(data_file,unpack=True))[6:]

    x_bin = bins_mean
    y_bin = profile_subsatboost[2]
    xerrlo_bin =  bins_low
    xerrhi_bin =  bins_up
    yerr_bin = profile_subsatboost[3]

    return(x_bin, y_bin, xerrlo_bin, xerrhi_bin, yerr_bin, lM_200m)
    

def Zhang_2024b_eRASS4_LX():

    cosmo = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)

    RES = np.transpose([
        [11.5, 12.0, 11.3, 11.8, 1.3, 0.6, 40, 5.4, 0.6, 39, 3.0, 0.9, 39, 2.4, 0.6, 39    ],
        [12.0, 12.5, 11.8, 12.3, 4.3, 0.7, 40, 1.9, 0.2, 40, 1.3, 0.2, 40, 6.3, 1.7, 39    ],
        [12.5, 13.0, 12.3, 12.8, 1.2, 0.2, 41, 6.5, 0.7, 40, 5.6, 0.7, 40, 8.5, 1.7, 39    ],
        [13.0, 13.5, 12.8, 13.3, 4.8, 1.1, 41, 2.8, 0.2, 41, 2.7, 0.2, 41, 1.4, 0.2, 40    ],
        [13.5, 14.0, 13.3, 13.7, 1.4, 0.1, 42, 1.1, 0.04, 42, 1.0, 0.04, 42, 3.2, 0.7, 40  ] ])
    M200min, M200max, M500min, M500max, LXtot_val, LXtot_err, LXtot_10e, LXmask_val, LXmask_err, LXmask_10e, LXCGM_val, LXCGM_err, LXCGM_10e, LXXPS_val, LXXPS_err, LXXPS_10e =RES
    zsample = np.array([0.08, 0.13, 0.16, 0.2, 0.2])
    M500c_mean = (M500max+M500min)*0.5
    #    plt.errorbar( 10**M500c_mean*cosmo.efunc(zsample), LXCGM_val*10**LXCGM_10e/cosmo.efunc(zsample),
    #            yerr = LXCGM_err*10**LXCGM_10e/cosmo.efunc(zsample),
    #            xerr = [ 10**M500c_mean*cosmo.efunc(zsample)-10**M500min*cosmo.efunc(zsample),10**M500max*cosmo.efunc(zsample)-10**M500c_mean*cosmo.efunc(zsample)],
    #            lw=3, ls='', color='purple', label='Zhang et al. 2024' )

    # return- M500, LXCGM, xerrlow,xerrhi, yerr 

    #print("Zhang LX errors= ", 10**M500c_mean*cosmo.efunc(zsample),10**M500min*cosmo.efunc(zsample),10**M500max*cosmo.efunc(zsample), LXCGM_err*10**LXCGM_10e/cosmo.efunc(zsample)) 
    return(10**M500c_mean*cosmo.efunc(zsample), LXCGM_val*10**LXCGM_10e/cosmo.efunc(zsample),10**M500c_mean*cosmo.efunc(zsample)-10**M500min*cosmo.efunc(zsample),10**M500max*cosmo.efunc(zsample)-10**M500c_mean*cosmo.efunc(zsample),LXCGM_err*10**LXCGM_10e/cosmo.efunc(zsample))
