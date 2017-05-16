from fabric.api import *



def pack():
    # build the package
    local('python setup.py bdist_wheel', capture=False)
    local('python setup.py sdist', capture=False)

def deploy():
    # update the pipy ham package
    run('twine upload dist/* -r testpypi')





