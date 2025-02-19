# DPMhalo

The Descriptive Parametric Model for gaseous halos version 1.0 updated on February 19, 2025

# Installation



via pip:

pip install dpmhalo

via git: 

git clone https://github.com/benopp99/DPMhalo.git

cd DPMhalo

python setup.py install

# Example Script



Upon installation, test that you can access a script to create halo radial profiles for a given M200c and redshift, via the script in your [site-packages directory]/dpmhalo/scripts/run_DPM_profiles.py.  

If you do not know your site-pacakges directory, you can obtain this script via:

https://github.com/benopp99/DPM/tree/main/dpmhalo/scripts/run_DPM_profiles.py

This script takes 2 or 3 arguments.  Run it once first with no arguments (python run_DPM_profiles.py) and it will explain the arguments of this script.  

For other scripts, go to:

https://github.com/benopp99/DPM/tree/main/dpmhalo/scripts/

There are a couple .csh shell scripts, and you need to do chmod a+x *.csh for them to become executables.  
