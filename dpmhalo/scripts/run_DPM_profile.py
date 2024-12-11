import numpy as np
import sys
import dpmhalo
from dpmhalo import myutils,darkmatter,modelprofiles

ModelDir = './'

if(len(sys.argv) < 3):
  print("Usage: run_DPM_profile.py 1) model name, 2) log(M200c), Optional: 3) redshift (default 0.0)\n model names: Model1, Model2, Model3, Model1ldsip, Model2ldisp, Model3ldisp, Model3disp\n")
  exit()

modelID = sys.argv[1]
lMhalo = float(sys.argv[2])
Mhalo = 10**lMhalo
if(len(sys.argv)>3):
  redshift = float(sys.argv[3])
else:
  redshift = 0.0

if(modelID.startswith("ClusterBased5") | modelID.startswith("Model1")):
  Pnorm12 = 4.4e+02
  alphatrP12 = 1.3
  alphahiP12 = 4.1
  alphaloP12 = 0.3
  c500P = 1.8
  alphatrPMvar = 0.0
  alphahiPMvar = 0.0 
  alphaloPMvar = 0.0
  betaP = 2/3.
  gammaP = 8/3.
  nenorm12 = 6.3e-04
  alphatrne12 = 1.9
  alphahine12 = 2.7
  alphalone12 = 1.0
  c500ne = 1.8
  alphatrneMvar = 0.0
  alphahineMvar = 0.0
  alphaloneMvar = 0.0
  betane = 0.0
  gammane = 2. 
  sigmalogne = 0.01
  Znorm12 = 0.2
  alphatrZ12 = 0.5
  alphahiZ12 = 0.7
  alphaloZ12 = 0.0
  c500Z = 1.8
  alphatrZMvar = 0.0
  alphahiZMvar = 0.0
  alphaloZMvar = 0.0
  betaZ = 0.0
  gammaZ = 0.0

if(modelID.startswith("ClusterScaled5") | modelID.startswith("Model2")):
  Pnorm12 = 1.24e+02
  alphatrP12 = 1.3
  alphahiP12 = 4.1
  alphaloP12 = 0.3
  c500P = 1.8
  alphatrPMvar = 0.0
  alphahiPMvar = 0.0 
  alphaloPMvar = 0.0
  betaP = 0.85
  gammaP = 8/3.
  nenorm12 = 5.24e-05
  alphatrne12 = 1.9
  alphahine12 = 2.7
  alphalone12 = 1.0
  c500ne = 1.8
  alphatrneMvar = 0.0
  alphahineMvar = 0.0
  alphaloneMvar = 0.0
  betane = 0.36
  gammane = 2. 
  sigmalogne = 0.01
  Znorm12 = 0.2
  alphatrZ12 = 0.5
  alphahiZ12 = 0.7
  alphaloZ12 = 0.0
  c500Z = 1.8
  alphatrZMvar = 0.0
  alphahiZMvar = 0.0
  alphaloZMvar = 0.0
  betaZ = 0.0
  gammaZ = 0.0

if(modelID.startswith("ClusterGroupScaled5") | modelID.startswith("Model3")):
  Pnorm12 = 7.6e+01
  alphatrP12 = 0.2
  alphahiP12 = 2.0
  alphaloP12 = -0.6
  c500P = 1.8
  alphatrPMvar = 0.3667
  alphahiPMvar = 0.7
  alphaloPMvar = 0.3
  betaP = 0.92
  gammaP = 8/3. 
  nenorm12 = 5.24e-05
  alphatrne12 = 0.45
  alphahine12 = 0.50
  alphalone12 = 0.40
  c500ne = 1.8
  alphatrneMvar = 0.483
  alphahineMvar = 0.733
  alphaloneMvar = 0.2
  betane = 0.36
  gammane = 2. 
  sigmalogne = 0.01
  Znorm12 = 0.2
  alphatrZ12 = 0.5
  alphahiZ12 = 0.7
  alphaloZ12 = 0.0
  c500Z = 1.8
  alphatrZMvar = 0.0
  alphahiZMvar = 0.0
  alphaloZMvar = 0.0
  betaZ = 0.0
  gammaZ = 0.0

if(modelID.endswith("disp")):
  sigmalogne = 0.30
if(modelID.endswith("ldisp")):
  sigmalogne = 0.15

ModelMFlexGNFWParams = dpmhalo.ModelMFlexGNFW(Mhalo,redshift,Pnorm12,alphatrP12,alphahiP12,alphaloP12,c500P,alphatrPMvar,alphahiPMvar,alphaloPMvar,betaP,gammaP,nenorm12,alphatrne12,alphahine12,alphalone12,c500ne,alphatrneMvar,alphahineMvar,alphaloneMvar,betane,gammane,sigmalogne,Znorm12,alphatrZ12,alphahiZ12,alphaloZ12,c500Z,alphatrZMvar,alphahiZMvar,alphaloZMvar,betaZ,gammaZ)

fout = open("%s/ModelMFlexGNFW%s.M%5.2f.z%4.2f.dat"%(ModelDir,modelID,lMhalo,redshift),"w")
fout.write("#lM200c z R/R200c P n_e[cm^-3] T[K] Z DM[pc/cm^3] tau_SZ y_SZ SB_SoftXray[erg/s/kpc^2] N_OVI[cm^-2] N_OVII[cm^-2] N_OVIII[cm^-2] sigma_ne rho_DM[g/cm^3]\n")

R200c_kpc = myutils.R200c_from_lM200c(np.log10(Mhalo),redshift)

for i in range(0,24,1): 
  R = 10**(-2.0+i*0.1)
  ne = ModelMFlexGNFWParams.calcne(R,Mhalo,redshift)
  x_kpc = R*R200c_kpc
  maxR_kpc = 3.0*R200c_kpc

  if(x_kpc<10.0/(1+redshift)): continue
  
  DM = ModelMFlexGNFWParams.AbelTransform("FRB",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  ySZ = ModelMFlexGNFWParams.AbelTransform("tSZ",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  tauSZ = ModelMFlexGNFWParams.AbelTransform("kSZ",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  SoftXray = ModelMFlexGNFWParams.AbelTransform("SoftXray",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  NOVI = ModelMFlexGNFWParams.AbelTransform("NOVI",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  NOVII = ModelMFlexGNFWParams.AbelTransform("NOVII",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  NOVIII = ModelMFlexGNFWParams.AbelTransform("NOVIII",x_kpc,maxR_kpc,R200c_kpc,Mhalo,redshift)
  
  P = ModelMFlexGNFWParams.calcP(R,Mhalo,redshift)
  Z = ModelMFlexGNFWParams.calcZ(R,Mhalo,redshift)
  Pe = P
  T = Pe/ne

  rho_DM = darkmatter.return_DMrho_for_R(R,np.log10(Mhalo),redshift)

  print("%5.2f %5.2f %5.3e %5.3e %5.3e %5.3e %5.3e %5.1f %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %4.2f %5.3e"%(np.log10(Mhalo),redshift,R,P,ne,T,Z,DM.value,tauSZ.value,ySZ.value,SoftXray.value,NOVI.value,NOVII.value,NOVIII.value,sigmalogne,rho_DM))
  fout.write("%5.2f %5.2f %5.3e %5.3e %5.3e %5.3e %5.3e %5.1f %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %4.2f %5.3e\n"%(np.log10(Mhalo),redshift,R,P,ne,T,Z,DM.value,tauSZ.value,ySZ.value,SoftXray.value,NOVI.value,NOVII.value,NOVIII.value,sigmalogne,rho_DM))

fout.close()
