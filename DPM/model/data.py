import numpy as np
import myconstants as myc 
from astropy import units as u
from astropy import constants as const
import smhm 
import myutils
from scipy.interpolate import RegularGridInterpolator

BaseDir = './'
    
def Arnaud_2010_Pe_fit_z0(lM200c,frac_lR200c):

    lM500c = myutils.lM200c_to_lM500c(lM200c)
    
    ###frac_R500c = myutils.R200_to_R500(10**frac_lR200c)
    frac_R500c = myutils.fR200c_to_fR500c(10**frac_lR200c, 10**lM200c)

    x = frac_R500c
    
    h_70 = 1.0
    Mass_exponent = (2/3.+(0.12+0.10)-(0.12+0.10)*(x/0.5)**3/(1+(x/0.5)**3))
    print("Mass_exponent, lM200c, frac_lR200c, x= ", Mass_exponent, lM200c, frac_lR200c, x)

    P_0 = 8.403*h_70**(-3/2.)
    c_500 = 1.177
    gamma = 0.3081
    alpha = 1.0510
    beta = 5.4905
 
    P_e = 1.65e-03*(10**lM500c/(3e+14/h_70))**Mass_exponent * P_0*h_70**2/((c_500*x)**gamma*(1+(c_500*x)**alpha)**((beta-gamma)/alpha))*myc.K_per_keV

    n_over_ne = 1/(myc.ne_to_nH*myc.XH*myc.mu)
    P = P_e*n_over_ne
    
    return(P_e)

def Ghirardini_2019_ne_fit(frac_lR200c):

    lM200c_Girardini2019 = 14.9 # XXXX CHECK
    ###frac_R500c = myutils.R200_to_R500(10**frac_lR200c)
    frac_R500c = myutils.fR200c_to_fR500c(10**frac_lR200c, 10**lM200c_Girardini2019)
    
    x = frac_R500c

    gamma = 3.
    n_0 = np.exp(-4.4) #10**-4.4
    r_core = np.exp(-3.0) #10**-3.0
    r_scale = np.exp(-0.29) #10**-0.29
    alpha = 0.89
    beta = 0.43
    epsilon = 2.86
    
    n_e_squared = n_0**2 * (x/r_core)**-alpha/(1.+x**2/r_core**2)**(3.*beta-alpha/2.)/(1.+x**gamma/r_scale**gamma)**(epsilon/gamma)
    n_e = np.sqrt(n_e_squared)
    print("n_e = ", n_e)
    print("x = ", x)
    return(n_e)

def Ghirardini_2019_Pe_fit(lM200c,frac_lR200c,redshift,norm_P500=False):

    lM500c = myutils.lM200c_to_lM500c(lM200c)
    ###frac_R500c = myutils.R200_to_R500(10**frac_lR200c)
    frac_R500c = myutils.fR200c_to_fR500c(10**frac_lR200c, 10**lM200c)

    h_70 = 1.0
    E_z = myutils.E(redshift) # XXXXX This is redshift = 0.  

    P_500 = myc.K_per_keV * 3.426e-03 * (10**lM500c/(h_70**-1*10**15.))**(2/3.)*E_z**(8/3.)*(myc.f_b/0.16)*(myc.mu/0.6)*(myc.mu_e/1.14)
    
    x = frac_R500c

    P_0 = 5.68
    c_500 = 1.49
    gamma = 0.43
    alpha = 1.33
    beta = 4.40

    Pe_over_P500 = P_0/(c_500*x)**gamma/(1+(c_500*x)**alpha)**((beta-gamma)/alpha)
    if(norm_P500):
        return(Pe_over_P500)
    else:
        return(Pe_over_P500*P_500)


def Sun_2011_P(cluster=False, return_P500=False):

    if(cluster):
        lRfrac, lP500_med, lP500_lo, lP500_hi = np.loadtxt("%s/Data/Sun2011_normalizedpressure.groups.dat"%BaseDir,usecols=(0,1,2,3),unpack=True)
        M500c_Sun2011 = 5e+14 #XXXX GUESS-FIX
    else:
        lRfrac, lP500_med, lP500_lo, lP500_hi = np.loadtxt("%s/Data/Sun2011_normalizedpressure.clusters.dat"%BaseDir,usecols=(0,1,2,3),unpack=True)
        M500c_Sun2011 = 7e+13
        
    #lRfrac_200 = np.log10(myutils.R500_to_R200(10**lRfrac))
    lRfrac_200 = np.log10(myutils.fR500c_to_fR200c(10**lRfrac,M500c_Sun2011))

    redshift_Sun2011 = 0.033
    
    lM200c_Sun2011 = myutils.lM500c_to_lM200c(np.log10(M500c_Sun2011))
    
    P500_Sun2011 = myutils.convert_Pe_to_Pe500(1.0,lM200c_Sun2011,0.0,return_only_P500=True)
    
    P_med = 10**lP500_med*P500_Sun2011
    P_lo = 10**lP500_lo*P500_Sun2011
    P_hi = 10**lP500_hi*P500_Sun2011

    if(return_P500):
        return(lRfrac_200,10**lP500_med,10**lP500_lo,10**lP500_hi)
    else:
        return(lRfrac_200,P_med,P_lo,P_hi)
    
def Lovisari_2015_ne_Groups():
    M500c_Lovisari_2015 = 10**13.55 # Need to check.  

    R500_scale_data, ne_med, ne_lo, ne_hi = np.loadtxt("%s/Data/Lovisari2015.ne_groups.dat"%BaseDir,usecols=(0,1,2,3),unpack=True)

    R200_scale_data = myutils.fR500c_to_fR200c(R500_scale_data,M500c_Lovisari_2015)
    
    return(R200_scale_data,ne_med,ne_lo,ne_hi)

def Lovisari_2015_M500_LX():
    h_70 = 1.0
    lM500_L15 = np.linspace(13.3,15.0,10)
    M500_L15 = 10**lM500_L15
    a_L15 = 1.39
    b_L15 = -0.12
    C1_L15 = 1e+43*h_70**-2.
    C2_L15 = 5e+13*h_70**-1.
    lLX_L15 = a_L15*np.log10(M500_L15/C2_L15)+b_L15+np.log10(C1_L15)

    return(lM500_L15, lLX_L15)


def McDonald_2017_ne_MH145_z0(lM200c,lRfrac,R200c_kpc,redshift):

    ###R500c_kpc = myutils.R200_to_R500(R200c_kpc)
    R500c_kpc = myutils.fR200c_to_fR500c(R200c_kpc, lM200c)

    M500c = 10**(myutils.lM200c_to_lM500c(lM200c))
    
    R_kpc = 10**lRfrac * R200c_kpc
    
    R_scale_500 = 10**lRfrac*(R500c_kpc/R200c_kpc) # scaled to R500, so multiply.
    R_scale_500_data, lrho_crit_500_data, err_500 = np.loadtxt("%s/Data/McDonald_2017.rhocrit_r500.dat"%BaseDir ,usecols=(0,1,2),unpack=True)
    lrho_crit_500 = np.interp(R_scale_500,R_scale_500_data,lrho_crit_500_data)

    rhocrit_cgs = 1.88e-29*myc.hubbleparam**2*(1+redshift)**3
    rhocrit = rhocrit_cgs/const.m_p.to('g').value

    nH_cgs = 10**lrho_crit_500*rhocrit*myc.XH

    ne_cgs = nH_cgs*myc.ne_to_nH

    R_kpc_return = np.linspace(10.,R200c_kpc,num=50)
    ne_cgs_return = np.interp(R_kpc_return,R_kpc,ne_cgs)

    return(R_kpc_return,ne_cgs_return)

def WuMcQuinn_2023_FRB(lMhalo):

    lmvirlo,lmvirhi,frvirlo,frvirhi,DM,DMerr = np.loadtxt("%s/Data/WuMcQuinn_2023.FRB.dat"%BaseDir,usecols=(0,1,2,3,4,5),unpack=True)
    lmvirmed = (lmvirlo+lmvirhi)/2.
    frvirmed = (frvirlo+frvirhi)/2.
    
    indexes = np.where((lMhalo>lmvirlo) & (lMhalo<lmvirhi))
    lmvirmed = lmvirmed[indexes]
    frvirmed = frvirmed[indexes]
    DM = DM[indexes]
    DMerr = DMerr[indexes]

    return(lmvirmed,frvirmed,DM,DMerr)

def Pratt_2021_tSZ():

    M500c_Pratt = 10**13.75
    R500_scale_data, ySZ_data, ySZ_err = np.loadtxt("%s/Data/Pratt_2021.groups_stack.dat"%BaseDir ,usecols=(0,1,2),unpack=True)
    ###R200_scale_data = myutils.R500_to_R200(R500_scale_data)
    R200_scale_data = myutils.fR500c_to_fR200c(R500_scale_data,M500c_Pratt)

    return(R200_scale_data,ySZ_data,ySZ_err)

def Bregman_2022_tSZ():

    # This uses median and assumes Rvir = R200.  
    R200_scale_data, ySZ_med, ySZ_err = np.loadtxt("%s/Data/Bregman_2022-updated.LStar_stack.dat"%BaseDir ,usecols=(1,3,5),delimiter=',',unpack=True)
    
    return(R200_scale_data,ySZ_med,ySZ_err)

def Schaan_2021_tSZ_CMASS():
    R_arcmin, y_CAP_data_ster, y_CAP_err_ster = np.loadtxt("%s/Data/data_schaan21/diskring_tsz_uniformweight_measured.txt"%BaseDir,usecols=(0,1,2),unpack=True)

    y_CAP_data_arcmin2 = y_CAP_data_ster*60**2*(180/np.pi)**2
    y_CAP_err_arcmin2 = y_CAP_err_ster*60**2*(180/np.pi)**2
    
    return(R_arcmin,y_CAP_data_arcmin2,y_CAP_err_arcmin2)

def Schaan_2021_kSZ_CMASS():
    R_arcmin_f090, TkSZ_data_f090, TkSZ_err_f090 = np.loadtxt("%s/Data/data_schaan21/f090/diskring_ksz_varweight_measured.txt"%BaseDir,usecols=(0,1,2),unpack=True)
    R_arcmin_f150, TkSZ_data_f150, TkSZ_err_f150 = np.loadtxt("%s/Data/data_schaan21/f150/diskring_ksz_varweight_measured.txt"%BaseDir,usecols=(0,1,2),unpack=True)

    tau_CAP_data_f090 = TkSZ_data_f090*60**2*(180/np.pi)**2/(2.76e+06*(313/2.998e+05))
    tau_CAP_err_f090 = TkSZ_err_f090*60**2*(180/np.pi)**2/(2.76e+06*(313/2.998e+05))

    tau_CAP_data_f150 = TkSZ_data_f150*60**2*(180/np.pi)**2/(2.76e+06*(313/2.998e+05))
    tau_CAP_err_f150 = TkSZ_err_f150*60**2*(180/np.pi)**2/(2.76e+06*(313/2.998e+05))

    return(R_arcmin_f090,tau_CAP_data_f090,tau_CAP_err_f090,R_arcmin_f150,tau_CAP_data_f150,tau_CAP_err_f150)


def Lovisari_2019_Z(combined=False):

    M500c_Lovisari_2019 = 10**13.55
    R500_lo_data, R500_hi_data, Z_relax_data, sig_relax_data, Z_dist_data, sig_dist_data = np.loadtxt("%s/Data/Lovisari_2019.Groups_Z.dat"%BaseDir, usecols=(0,1,2,3,4,5),unpack=True)

    Z_relax_data *= myc.Z_Solar_Asplund/myc.Z_Solar_Anders # normalize to Anders & Grevesse 1989 since in Asplund 2009
    Z_dist_data *= myc.Z_Solar_Asplund/myc.Z_Solar_Anders # normalize to Anders & Grevesse 1989 since in Asplund 2009
    
    ###R200_scale_data = myutils.R500_to_R200((R500_hi_data+R500_lo_data)/2.)
    R200_scale_data = myutils.fR500c_to_fR200c((R500_hi_data+R500_lo_data)/2.,M500c_Lovisari_2019)

    if(combined is True):
        Z_combined_data = (Z_relax_data+Z_dist_data)/2.
        sig_combined_data = (sig_relax_data+sig_dist_data)/2.
        return(R200_scale_data,Z_combined_data,sig_combined_data)
    else:
        #Return relaxed groups only
        return(R200_scale_data,Z_relax_data,sig_relax_data)

def Ghizzardi_2021_Z():
    M500c_Ghizzardi_2021 = 10**14.75 #XXX CHECK
    R500c_lo_data, R500c_hi_data, ZFe_med_data, ZFe_err_data = np.loadtxt("%s/Data/Ghizzardi_2021.Cluster_Z.dat"%BaseDir, usecols=(0,1,5,6),unpack=True)
    R500c_scale_data = (R500c_lo_data+R500c_hi_data)/2.
    ###R200_scale_data = myutils.R500_to_R200(R500_scale_data)
    R200_scale_data = myutils.fR500c_to_fR200c(R500c_scale_data,M500c_Ghizzardi_2021)

    return(R200_scale_data,ZFe_med_data,ZFe_err_data)

def XGAP_profiles(lM200c, return_ne=True,return_P=False,return_ZFe=False):

    lM500c_read, ngal, frac_lR500c_lo, frac_lR500c_hi, lPmed, lPlo, lPhi, lnemed, lnelo, lnehi, lTmed, lTlo, lThi = np.loadtxt("%s/Data/XGAP_median_profiles.dat"%BaseDir, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12), unpack=True)
    frac_lR500c = (frac_lR500c_lo+frac_lR500c_hi)/2.

    lPmed = lPmed + np.log10(myc.K_per_keV)
    lPlo = lPlo + np.log10(myc.K_per_keV)
    lPhi = lPhi + np.log10(myc.K_per_keV)

    lM500c = myutils.lM200c_to_lM500c(lM200c)
    if((lM500c>=13.2) & (lM500c<13.53333)):
        lM500clo = 13.2
        lM500chi = 13.53333
    if((lM500c>=13.5333) & (lM500c<13.86667)):
        lM500clo = 13.53333
        lM500chi = 13.86667
    if((lM500c>=13.86667) & (lM500c<=14.2)):
        lM500clo = 13.86667
        lM500chi = 14.2

    indexes_return = np.where((lM500c_read>=lM500clo) & (lM500c_read<lM500chi))

    ###frac_R200c = myutils.R500_to_R200(10**frac_lR500c)
    frac_R200c = myutils.fR500c_to_fR200c(10**frac_lR500c,10**lM500c)

    if(return_ZFe):
        lM500c_read, ngal, frac_R500c_lo, frac_R500c_hi, lZmed, lZlo, lZhi = np.loadtxt("%s/Data/XGAP_median_Fe_profiles.dat"%BaseDir, usecols=(0,1,2,3,4,5,6), unpack=True)
        frac_R500c = (frac_R500c_lo+frac_R500c_hi)/2.
        indexes_return = np.where((lM500c_read>=lM500clo) & (lM500c_read<lM500chi))

        ###frac_R200c = myutils.R500_to_R200(frac_R500c)
        frac_R200c = myutils.fR500c_to_fR200c(frac_R500c,10**lM500c)

        return(frac_R200c[indexes_return],lZmed[indexes_return],lZlo[indexes_return],lZhi[indexes_return])
    if(return_P):
        return(frac_R200c[indexes_return],lPmed[indexes_return],lPlo[indexes_return],lPhi[indexes_return])
    if(return_ne):
        return(frac_R200c[indexes_return],lnemed[indexes_return],lnelo[indexes_return],lnehi[indexes_return])

def Voit_logU_to_logP(logU,redshift):

    redshift_vec = [ 0.0 , 0.4 ]
    log_nph_vec = [ -5.97 , -5.42 ]
    log_nph  = np.interp(redshift, redshift_vec, log_nph_vec)

    nH = 10**log_nph/10.0**logU

    logU_vec = np.flip([-1.42, -1.67, -1.92, -2.17, -2.42, -2.67, -2.92, -3.17, -3.42, -3.67, -3.92, -4.17, -4.42, -4.67, -4.92, -5.17, -5.42 ])
    logT_vec = np.flip([ 4.342, 4.291, 4.291, 4.192, 4.146, 4.105, 4.071, 4.044, 4.022, 4.000, 3.977, 3.954, 3.932, 3.911, 3.888, 3.861, 3.827])
    log_T = np.interp(logU, logU_vec, logT_vec)

    T = 10**log_T
    ne = nH*myc.ne_to_nH
    
    Ptherm = T*ne

    return(np.log10(Ptherm))
    
def Voit_lognH_to_logP(lognH,redshift):

    redshift_vec = [ 0.0 , 0.4 ]
    log_nph_vec = [ -5.97 , -5.42 ]
    log_nph  = np.interp(redshift, redshift_vec, log_nph_vec)

    U = 10**log_nph/10**lognH
    logU = np.log10(U)
    
    logU_vec = np.flip([-1.42, -1.67, -1.92, -2.17, -2.42, -2.67, -2.92, -3.17, -3.42, -3.67, -3.92, -4.17, -4.42, -4.67, -4.92, -5.17, -5.42 ])
    logT_vec = np.flip([ 4.342, 4.291, 4.291, 4.192, 4.146, 4.105, 4.071, 4.044, 4.022, 4.000, 3.977, 3.954, 3.932, 3.911, 3.888, 3.861, 3.827])
    log_T = np.interp(logU, logU_vec, logT_vec)

    T = 10**log_T
    ne = 10**lognH*myc.ne_to_nH

    Ptherm = T*ne
    return(np.log10(Ptherm))

def Voit_LogU(sample=False):

    if(sample == "COS-Halos"):
        logU,logU_hi,logU_lo,Rproj,lMstar_pub,redshift,fracRvir = np.loadtxt("%s/Data/LogU_table.COS-Halos.csv"%BaseDir, usecols=(2,3,4,5,6,7,8), delimiter=',',unpack=True)
    else:
        if(sample == "no-COS-Halos"):
            logU,logU_hi,logU_lo,Rproj,lMStar_pub,redshift,fracRvir = np.loadtxt("%s/Data/LogU_table.no-COS-Halos.csv"%BaseDir, usecols=(2,3,4,5,6,7,8), delimiter=',',unpack=True)
        else:       
            logU,logU_hi,logU_lo,Rproj,lMStar_pub,redshift,fracRvir = np.loadtxt("%s/Data/LogU_table.csv"%BaseDir, usecols=(2,3,4,5,6,7,8), delimiter=',',unpack=True)

    lMStar_err = 0.2
    #lM200 = smhm.Behroozi2013_return_mhalo(lMStar_pub,0.0)
    #return(Rproj,logU,logU_hi,logU_lo,lM200,redshift)

    lM200 = lMStar_pub*0.
    lM200_errpos = lMStar_pub*0.
    lM200_errneg = lMStar_pub*0.
    for i in range(len(lMStar_pub)):
        lM200[i], lM200_errpos[i], lM200_errneg[i] = smhm.Behroozi2019_UNIVERSEMACHINE_return_mhalo(lMStar_pub[i],lMStar_err,redshift[i])

    return(Rproj,logU,logU_hi,logU_lo,lM200,lM200_errpos,lM200_errneg,redshift)

def Prochaska_2017_Densities():
    lMStar_pub, sSFR_pub, Rproj, lognH_lo, lognH, lognH_hi = np.loadtxt("%s/Data/Prochaska_2017.COS-Halos.dat"%BaseDir, usecols=(2,3,4,5,6,7), unpack=True)
    # We should read redshift, but it's weird... not really needed for lognH, unlike logU.   
    redshift = lMStar_pub*0 + 0.2
    lMStar_err = 0.2
    
    #lM200 = smhm.Behroozi2013_return_mhalo(lMStar_pub,0.0)
    #return(Rproj, lognH, lognH_lo, lognH_hi, lM200, redshift)

    lM200 = lMStar_pub*0.
    lM200_errpos = lMStar_pub*0.
    lM200_errneg = lMStar_pub*0.
    for i in range(len(lMStar_pub)):
        lM200[i], lM200_errpos[i], lM200_errneg[i] = smhm.Behroozi2019_UNIVERSEMACHINE_return_mhalo(lMStar_pub[i],lMStar_err,redshift[i])
        
    return(Rproj, lognH, lognH_lo, lognH_hi, lM200,lM200_errpos,lM200_errneg, redshift)

def Tchernyshyov_2022_NOVI(survey_return):

    types = ["|S12", "float", "float", "float", "|S2", "float", "float","|S1"]

    survey,redshift,R,lMStar, MType,R200c,lNOVI,OVI_flag = np.genfromtxt("%s/Data/Tchernyshyov_2022.Table3.tsv"%BaseDir,delimiter=';',dtype=types,usecols=(0,4,5,6,7,8,9,10),unpack=True)
    #for i in range(len(survey)):
    #    survey[i] = str(survey[i], "utf-8")
    #    #survey[i] = survey[i].decode("utf-8")
    #print("OVI_flag= ", OVI_flag)
    #print("survey= ", survey)

    lMStar_err = 0.2
    
    R_return = []
    redshift_return = []
    R200c_return = []
    lMStar_return = []
    lNOVI_return = []
    OVI_flag_int_return = []
    lM200c_smhm_return = []
    lM200c_errpos_smhm_return = []
    lM200c_errneg_smhm_return = []
    survey_str = [None]*len(lNOVI)
    OVI_flag_int = np.zeros(len(lNOVI))
    for i in range(len(survey)):
        survey_str[i] = survey[i].decode("utf-8")
        if(OVI_flag[i].decode("utf-8")=="<"):
            OVI_flag_int[i] = 1
        else:
            OVI_flag_int[i] = 0
        if survey_return in survey_str[i]:
            #print(redshift[i],R[i]/R200c[i],lMStar[i])
            lM200c_smhm, lM200c_errpos_smhm, lM200c_errneg_smhm = smhm.Behroozi2019_UNIVERSEMACHINE_return_mhalo(lMStar[i],lMStar_err,redshift[i])
            redshift_return.append(redshift[i])
            R_return.append(R[i])
            lMStar_return.append(lMStar[i])
            R200c_return.append(R200c[i])
            lNOVI_return.append(lNOVI[i])
            OVI_flag_int_return.append(OVI_flag_int[i])
            lM200c_smhm_return.append(lM200c_smhm)
            lM200c_errpos_smhm_return.append(lM200c_errpos_smhm)
            lM200c_errneg_smhm_return.append(lM200c_errneg_smhm)
            

    #return(np.asarray(R_return), np.asarray(R200c_return), np.asarray(lNOVI_return), np.asarray(OVI_flag_int_return), np.asarray(lMStar_return), np.asarray(redshift_return))
    return(np.asarray(R_return), np.asarray(R200c_return), np.asarray(lNOVI_return), np.asarray(OVI_flag_int_return), np.asarray(lM200c_smhm_return), np.asarray(lM200c_errpos_smhm_return), np.asarray(lM200c_errpos_smhm_return), np.asarray(redshift_return))
    



def Miller_2015_nH_fit_MW(lM200c,lR,R200c):

    R200c_kpc = R200c*1e+03
    r_kpc = 10**lR * R200c_kpc
    ne_over_nh = 1.16
    
    n0_rc_normed = 1.35e-02
    correction = 1./0.3
    #correction = 1.1e-04/1.54e-05 # Corrected to Faerman LMC e^-2 value. 
    beta= 0.50
    
    n_e = n0_rc_normed*correction/(r_kpc)**(3*beta)

    print("r_kpc = ", r_kpc)
    print("n_e = ", n_e)

    nH = n_e/ne_over_nh

    M_200_hot = myutils.calc_spherical_mass(nH,r_kpc,R200c_kpc)
    print("M_200_hot_Miller = ", M_200_hot)
    return(nH)

def Bregman_2018_nH_fit_MW(lM200c,lR,R200c):
    #n(r) = n_0 r_c^(3*beta)/r^(3*beta)
    #n_0 r_c^(3*beta) = 1.20e-02 (+2.13,-0.82)
    #beta = 0.56 (+0.10,-0.12)
    # Is this electron density?  XXX

    R200c_kpc = R200c*1e+03
    r_kpc = 10**lR * R200c_kpc
    ne_over_nh = 1.16
    
    n0_rc_normed = 1.2e-02
    correction = 1.1e-04/1.54e-05 # Corrected to Faerman LMC e^-2 value. 
    beta= 0.56
    
    n_e = n0_rc_normed*correction/(r_kpc)**(3*beta)

    print("r_kpc = ", r_kpc)
    print("n_e = ", n_e)

    nH = n_e/ne_over_nh

    M_200_hot = myutils.calc_spherical_mass(nH,r_kpc,R200c_kpc)
    print("M_200_Bregman2018 = ", M_200_hot)
    return(nH)

