from fabric.api import *


def clean():
    # remove old package
    local('rm -r dist/', capture=False)
    local('rm -r build/', capture=False)

def test():
    # run unit test
    local('python -m unittest discover tests/', capture=False)
    local('python tests/functional_test.py ', capture=False)

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





