from fabric.api import *


def clean():
    # remove old package
    local('rm -r dist/', capture=False)
    local('rm -r build/', capture=False)

def pack():
# build the package
    local('python setup.py bdist_wheel', capture=False)
    local('python setup.py sdist', capture=False)

def deploy_test():
    # update the pipy ham package
    local('twine upload dist/* -r testpypi')





