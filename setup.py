#from distutils.core import setup, Extension
#import os
from setuptools import setup


setup(name='DPM',
      version='0.1',
      description='The Descriptive Parametric Model',
      url='https://github.com/benopp99/DPM',
      author='Benjamin D. Oppenheimer',
      author_email='beop5934@colorado.edu',
      license='MIT',
      packages=['DPM'],
      #package_dir={'DPM':'dpm'},
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

