======================================================================
trustregion: Trust-region subproblem solver for nonconvex optimization
======================================================================

.. image::  https://travis-ci.org/lindonroberts/trust-region.svg?branch=master
   :target: https://travis-ci.org/lindonroberts/trust-region
   :alt: Build Status

.. image::  https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: GNU GPL v3 License

General information here

Documentation
-------------


Citation
--------

Requirements
------------
:code:`trustregion` requires the following software to be installed:

* Python 2.7 or Python 3 (http://www.python.org/)

Additionally, the following python packages should be installed (these will be installed automatically if using *pip*, see `Installation using pip`_):

* NumPy 1.11 or higher (http://www.numpy.org/)

Installation using pip
----------------------
For easy installation, use `pip <http://www.pip-installer.org/>`_ as root:

 .. code-block:: bash

    $ [sudo] pip install numpy
    $ [sudo] pip install trustregion

Note that NumPy should be installed before :code:`trustregion`, as it is used to compile the Fortran modules.

If you do not have root privileges or you want to install :code:`trustregion` for your private use, you can use:

 .. code-block:: bash

    $ pip install --user numpy
    $ pip install --user trustregion

which will install :code:`trustregion` in your home directory.

Note that if an older install of :code:`trustregion` is present on your system you can use:

 .. code-block:: bash

    $ [sudo] pip install --upgrade trustregion

to upgrade :code:`trustregion` to the latest version.

Manual installation
-------------------
Alternatively, you can download the source code from `Github <https://github.com/lindonroberts/trust-region>`_ and unpack as follows:

 .. code-block:: bash

    $ git clone https://github.com/lindonroberts/trust-region
    $ cd trust-region

To upgrade :code:`trustregion` to the latest version, navigate to the top-level directory (i.e. the one containing :code:`setup.py`) and rerun the installation using :code:`pip`, as above:

 .. code-block:: bash

    $ git pull
    $ [sudo] pip install .  # with admin privileges

Testing
-------
If you installed :code:`trustregion` manually, you can test your installation by running:

 .. code-block:: bash

    $ python setup.py test

Alternatively, the documentation provides some simple examples of how to run :code:`trustregion`.

Examples
--------
Examples of how to run :code:`trustregion` are given in the `documentation <https://lindonroberts.github.io/trust-region/>`_, and the `examples <https://github.com/lindonroberts/trust-region/tree/master/examples>`_ directory in Github.

Uninstallation
--------------
If :code:`trustregion` was installed using *pip* you can uninstall as follows:

 .. code-block:: bash

    $ [sudo] pip uninstall trustregion

If :code:`trustregion` was installed manually you have to remove the installed files by hand (located in your python site-packages directory).

Bugs
----
Please report any bugs using GitHub's issue tracker.

License
-------
This algorithm is released under the GNU GPL license.
