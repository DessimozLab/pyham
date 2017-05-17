How to install pyham ?
====================

.. note:: We encourage pyham users to work on virtual environment. If you don't understand the last sentences you should consider reading the following short guide to python virtual environment: http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/

Installing via pip
##################

pyham can be install via pip using the following command:


    - **First upgrade pip to the latest version**:
        .. code-block:: bash

                python -m pip install --upgrade pip

    - **Then install the pyham package**:

        .. code-block:: bash

                pip install pyham

        in case you don't have root access (and/or permissions for writing in system directories) you can use the --user flag to install in the local user package  :

        .. code-block:: bash

                pip install --user pyham.

Source code
###########

The lastest source code of pyham version |release| can be download here: **www.lab.dessimoz/ham/latest.tar.gz**.


Install PyQt4 on Mac (only required for treeProfile)
####################################################


**Thanks to Dr. D. Dylus for writing this section.**


Installing PyQt4 on mac with pyenv adapted from: https://wiki.auckland.ac.nz/display/CERES/Installing+Qt+and+PyQt+libraries+for+ete3+-+tree+visualization+library

Below are my current system specs:
Mac OsX 10.9.5
Python 3.5.1

To get it working, you need to install these libraries in a particular order:

 - Install Qt
 - Install SIP
 - Install PyQt
 - Install ete3


--------------------------------------


Install Qt
----------

**Don't download the latest version of Qt**, which is 5.5 (this is the latest version when I installed it). The libraries (for example, ete3) built on previous versions of Qt will have issues with this latest version of Qt, as it requires you to include QtDBus.framework in your app, or otherwise it crashes. So, the libraries built on previous versions haven't included the framework and crashes with a false error:
This application failed to start because it could not find or load the Qt platform plugin "cocoa".
You can see an issue raised for this in Qt forums.

Download the Qt version 5.4.2 from http://download.qt.io/archive/qt/5.4/5.4.2/

Just double click the file, and it will open a setup box. Just follow it. You will be required to create a Qt account, if you don't have one. Qt is important since qmake is essential in order to compile pyqt4. You can choose the default installation folder or provide a location where you normally installed all your packages (in my case /Users/name/.pyenv/versions/3.5.1/envs/notebook3/lib/python3.5/site-packages).


--------------------------------------

Install sip
-----------
Download the current version of SIP from

https://www.riverbankcomputing.com/software/sip/download

.. code-block:: bash

                python configure.py -d /Users/name/.pyenv/versions/3.5.1/envs/notebook3/lib/python3.5/site-packages

Finally, make it and install it with the following commands:


.. code-block:: bash

                make
                sudo make install


--------------------------------------

Install PyQt4
-------------

Download the mac version from https://www.riverbankcomputing.com/software/pyqt/download

.. code-block:: bash

                python configure.py -q /Users/name/.pyenv/versions/3.5.1/envs/notebook3/lib/python3.5/site-packages/5.4/clang_64/bin/qmake -d /Users/name/.pyenv/versions/3.5.1/envs/notebook3/lib/python3.5/site-packages

if you have a problem with qmake on this stage do the following:

.. code-block:: bash

                cd /Applications/Xcode.app/Contents/Developer/usr/bin/
                sudo ln -s xcodebuild xcrun

then repeat the configure.py

.. code-block:: bash

                make
                sudo make install

--------------------------------------

Some additional links on the topics that can help:
 - Tutorial on installing PyQt and its dependency SIP: http://movingthelamppost.com/blog/html/2013/07/12/installing_pyqt____because_it_s_too_good_for_pip_or_easy_install_.html
 - Tutorial on installing PySide and PyQt on Windows, Mac and Linux: http://pythoncentral.io/install-pyside-pyqt-on-windows-mac-linux/
