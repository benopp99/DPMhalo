# DPMhalo

The Descriptive Parametric Model for gaseous halos version 0.1 created on August 28, 2024

# Brief Description

This is a project to create atmospheric radial profiles for gaseous halos as a function of halo mass (M) and redshift (z).  This model only creates 1-D spherical profiles and is meant to model the volume-filling component of the circumgalactic medium (CGM).  Embedded, cloud-like structures are not directly modeled.  A paper will be released, Oppenheimer, B. D.; Voit, G. M. et al..  In the meantime, you are welcome to download version 0.X of this code.  We encourage users to use these profiles in their research.  

# Installation



via pip:

pip install git+https://github.com/benopp99/DPMhalo.git#egg=dpmhalo



via git: 

git clone https://github.com/benopp99/DPMhalo.git

cd dpmhalo

python setup.py install

# Example Script



Upon installation, test that you can access a script to create halo radial profiles for a given M200c and redshift, via the script in your [site-packages directory]/dpm/scripts/run_DPM_profiles.py.  

If you do not know your site-pacakges directory, you can obtain this script via:

https://github.com/benopp99/DPM/tree/main/dpm/scripts/run_DPM_profiles.py

This script takes 2 or 3 arguments.  Run it once first with no arguments (python run_DPM_profiles.py) and it will explain the arguments of this script.  



