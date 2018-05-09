from fabric.api import *


# TO UPDATE PIP PACKAGE
# 1 -> update __version__ variable in pyham/__init__.py
# 2 -> fab test clean pack deploy


def clean():
    # remove old package
    local('rm -rf dist/', capture=False)
    local('rm -rf build/', capture=False)

def test():
    # run unit test
    local('python -m unittest discover tests/', capture=False)
    local('python tests/functional_test.py', capture=False)

def pack():
    # build the package
    local('python setup.py bdist_wheel', capture=False)
    local('python setup.py sdist', capture=False)

def deploy_test():
    # update the pipy pyham package
    local('twine upload dist/* -r testpypi')

def deploy():
    # update the pipy pyham package
    local('twine upload dist/* -r pypi')





