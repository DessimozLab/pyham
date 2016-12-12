from setuptools import setup
import os

name = 'ham'

setup(
    name=name,
    version='0.0.0',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='A tool to analyse Hierarchical Orthologous Groups (HOGs)',
    install_requires=['lxml'],
)
