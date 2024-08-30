import numpy as np
import h5py
import smhm 
import myutils

ModelDir = '.'

def read_COSHalos():

    COSHalos_path = ModelDir + "/Data/absorbers/COS-Halos/COS-Halos.hdf5"
    COSHalos_in = h5py.File(COSHalos_path)
    
    #EW_Lya_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['EW_Lya'])
    #EW_Lya_tag_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['EW_Lya_tag'], dtype= np.int )
    #EW_Lya_err_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['EW_Lya_err'])
    #EW_SiIII_COSHalos = np.asarray(COSHalos_in['Pairs']['SiIII']['EW_SiIII'])
    #EW_SiIII_tag_COSHalos = np.asarray(COSHalos_in['Pairs']['SiIII']['EW_SiIII_tag'], dtype= np.int )
    #EW_SiIII_err_COSHalos = np.asarray(COSHalos_in['Pairs']['SiIII']['EW_SiIII_err'])
    #print EW_SiIII_COSHalos
    #NHI_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['N_HI'])
    #NHI_tag_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['N_HI_tag'], dtype= np.int )
    #NHI_err_COSHalos = np.asarray(COSHalos_in['Pairs']['HI']['N_HI_err'])
    NOVI_COSHalos = np.asarray(COSHalos_in['Pairs']['OVI']['N_OVI'])
    NOVI_tag_COSHalos = np.asarray(COSHalos_in['Pairs']['OVI']['N_OVI_tag'])
    NOVI_err_COSHalos = np.asarray(COSHalos_in['Pairs']['OVI']['N_OVI_err'])
    bgal_COSHalos = np.asarray(COSHalos_in['Pairs']['bgal'])
    zgal_COSHalos = np.asarray(COSHalos_in['Pairs']['zgal'])
    lmstar_COSHalos = np.asarray(COSHalos_in['Pairs']['lmstar'])
    lssfr_COSHalos = np.asarray(COSHalos_in['Pairs']['lssfr'])

    lmstar_err = 0.3    
    redshift = 0.2 # XXXX NEED ACTUAL REDSHIFT XXX
    
    #lM200_COSHalos = smhm.Behroozi2013_return_mhalo(lmstar_COSHalos,0.0)

    lM200_COSHalos = lM200_errpos_COSHalos = lM200_errneg_COSHalos = lmstar_COSHalos*0.
    for i in range(len(lmstar_COSHalos)):
        lM200_COSHalos[i], lM200_errpos_COSHalos[i], lM200_errneg_COSHalos[i] = smhm.Behroozi2019_UNIVERSEMACHINE_return_mhalo(lmstar_COSHalos[i],lmstar_err,redshift)

    #print("lM200_COSHalos= ", lM200_COSHalos)
    #print("bgal_COSHalos= ", bgal_COSHalos)
    #print("NOVI_COSHalos= ", NOVI_COSHalos)

    return(bgal_COSHalos,NOVI_COSHalos,NOVI_tag_COSHalos,NOVI_err_COSHalos,lM200_COSHalos,lM200_errpos_COSHalos,lM200_errneg_COSHalos,zgal_COSHalos)


#read_COSHalos()
