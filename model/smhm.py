import numpy as np
import myutils

#This is the replacement for the Behroozi2013 code below and is includes returning error ranges.  
def Behroozi2019_UNIVERSEMACHINE_return_mhalo(lMStar,err_lMStar,redshift,galaxy_type_tag=''):

    redshift_tag = round(redshift*10)
    
    file_name = 'Data/Behroozi_UNIVERSEMACHINE/SMHM_med_cen%s.z0.%d.txt'%(galaxy_type_tag,redshift_tag)
    
    lmvir, lms, lms_errplus, lms_errneg = np.loadtxt(file_name,usecols=(0,1,3,4),unpack=True)

    lm200c = myutils.lMvir_to_lM200c(lmvir,redshift)
    
    # Sometimes the error return is negative and unrealistic, so we zero it.  
    lms_errplus = np.where(lms_errplus>0,lms_errplus,0)
    lms_errneg = np.where(lms_errneg>0,lms_errneg,0)

    lms_max = lms+lms_errplus
    lms_min = lms-lms_errneg

    lMHalo_med = np.interp(lMStar, lms, lm200c)

    # This adds errors linearly, not in quadarture, because we are making conservative (wide) errors.  
    lMHalo_max = np.interp(lMStar+err_lMStar, lms_min, lm200c)
    lMHalo_min = np.interp(lMStar-err_lMStar, lms_max, lm200c)

    #print("lMHalo_med, lMHalo_max-lMHalo_med, lMHalo_med-lMHalo_min= ", lMHalo_med, lMHalo_max-lMHalo_med, lMHalo_med-lMHalo_min)

    return(lMHalo_med, lMHalo_max-lMHalo_med, lMHalo_med-lMHalo_min)
    
#This is a very rudimentary code and has been replaced.  
def Behroozi2013_return_mhalo(lmstar_in,z):

    logM_step = 0.05
    logM_lo = 8.0
    logM_hi = 15.0
    logM = np.linspace(logM_lo,logM_hi,int((logM_hi-logM_lo)/logM_step)+1)
    
    mh = 10**logM
    a = 1./(1+z)

    eps0 = -1.785
    epsa = -0.074
    epsz = -0.048
    epsa2 = -0.179
    ml0 = 11.539
    mla = -1.751 
    mlz = -0.329
    alpha0 = -1.474
    alphaa = 1.339
    delta0 = 3.529
    deltaa = 4.152
    deltaz = 1.122
    gamma0 = 0.395
    gammaa = 0.766
    gammaz = 0.435
    nu = np.exp(-4*a**2)
    ml = 10**(ml0+(mla*(a-1)+mlz*z)*nu)
    eps = 10**(eps0 + (epsa*(a-1)+epsz*z)*nu + epsa2*(a-1))
    alpha = alpha0 + alphaa*(a-1)*nu
    delta = delta0 + (deltaa*(a-1)+deltaz*z)*nu
    gammag = gamma0 + (gammaa*(a-1)+gammaz*z)*nu

    x = np.log10(mh/ml)
    f =  -np.log10(10**(alpha*x)+1) + delta * np.log10(1+np.exp(x))**gammag / (1+np.exp(10**-x))
    f0 = -np.log10(10**(alpha*0)+1) + delta * np.log10(1+np.exp(0))**gammag / (1+np.exp(10**-0))
    mstar = 10**(np.log10(eps*ml) + f - f0)


    lmhalo_out = np.zeros(len(lmstar_in))
    for j in range(len(lmstar_in)):
        lmhalo_out[j] = np.interp(lmstar_in[j],np.log10(mstar),logM)
        #for i in range(len(logM)):
            #if(np.log10(mstar[i])>=lmstar_in[j]):
            #    lmhalo_out[j] = logM[i]
            #    lmhalo_out[j] = myutils.lMvir_to_lM200(lmhalo_out[j])
            #    #print("logM= ", logM[i],np.log10(mstar[i]))
            #    break

    return(lmhalo_out)
        