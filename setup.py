#from distutils.core import setup, Extension
#import os
from setuptools import find_packages, setup
from setuptools.extension import Extension


setup(name='DPMhalo',
      version='0.2',
      description='The Descriptive Parametric Model for gaseous halos',
      url='https://github.com/benopp99/DPMhalo',
      author='Benjamin D. Oppenheimer',
      author_email='beop5934@colorado.edu',
      license='MIT',
      #packages=['dpmhalo'],
      packages=find_packages(),
      include_package_data=True,
      #package_dir={'DPMhalo':'dpmhalo'},
      zip_safe=False,
      install_requires=['numpy',
                        'astropy',
                        'scipy',
                        'h5py',
                        'colossus',
                        'trident',
                        'pyxsim'
      ]
)


