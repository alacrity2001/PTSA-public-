from setuptools import setup
from os import path
from io import open
import setuptools

#this_directory = path.abspath(path.dirname(__file__))

#with open(path.join(this_directory, 'requirements.txt'),
#          encoding='utf-8') as f:
#    requirements = f.read().splitlines()

setup(name='ptsa',
      version='0.1',
      description='Python toolkil for timeseries anomoly detection',
      url='https://github.com/alacrity2001/SensorCleaning/tree/master/ptsa',
      author='Yinchen Wu',
      author_email='yinchen@uchicago.edu',
      license='MIT',
      packages=setuptools.find_packages(),
      #install_requires=requirements,
      zip_safe=False)