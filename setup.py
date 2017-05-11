from setuptools import setup, find_packages
import os

name = 'HAM'

__version__ = "Undefined"
for line in open('{}/__init__.py'.format(name.lower())):
    if line.startswith('__version__'):
        exec(line.strip())

requirements = ['ete3','six','scipy']

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=name,

    version=__version__,

    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',

    description='A tool to analyse Hierarchical Orthologous Groups (HOGs)',
    long_description=long_description,

    license='MIT',

    classifiers=[
         'Development Status :: 3 - Alpha',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'License :: OSI Approved :: MIT License',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Programming Language :: Python :: 3.4',
         'Programming Language :: Python :: 3.5',
         'Programming Language :: Python :: 3.6',
         ],

    packages=find_packages(),
    install_requires=requirements,
    extras_require={
        'test': ['sphinx','noise'],
    }
)
