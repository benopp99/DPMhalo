from colossus.cosmology import cosmology
cosmo = cosmology.setCosmology('planck15')
#my_cosmo = {'flat': True, 'H0': 67, 'Om0': 0.31, 'Ob0': 0.0490, 'sigma8': 0.8102, 'ns': 0.9665}
#cosmo = cosmology.setCosmology('my_cosmo', **my_cosmo)
##Cosmology "planck18" 
###    flat = True, Om0 = 0.3111, Ode0 = 0.6888, Ob0 = 0.0490, H0 = 67.66, sigma8 = 0.8102, ns = 0.9665
###    de_model = lambda, relspecies = True, Tcmb0 = 2.7255, Neff = 3.0460, powerlaw = False

#print("cosmo= ", cosmo)

hubbleparam = cosmo.H0/100. #0.7
HubbleParam = cosmo.H0  #100.*hubbleparam
OmegaM = cosmo.Om0 #0.3
OmegaL = 1-cosmo.Om0  #1.-OmegaM
Omegab = cosmo.Ob0 #0.048

ne_to_nH = 1.16
XH = 0.75
mu = 0.59
mu_e = 1./XH/ne_to_nH
f_b = Omegab/OmegaM
cm_per_kpc = 3.086e+21
km_per_Mpc = 3.086e+19
keV_per_K = 8.6e-08
K_per_keV = 1.16e+07

M_solar_g = 1.989e+33 

nO_per_nH_Zsolar = 0.00048977881 # 10**-3.31
Z_Solar_Asplund = 0.0134
Z_Solar_Anders = 0.0201
