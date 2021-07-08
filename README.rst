===========================================
trustregion: Trust-region subproblem solver
===========================================

.. image::  https://travis-ci.org/lindonroberts/trust-region.svg?branch=master
   :target: https://travis-ci.com/lindonroberts/trust-region
   :alt: Build Status

.. image::  https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: GNU GPL v3 License

.. image:: https://img.shields.io/pypi/v/trustregion.svg
   :target: https://pypi.python.org/pypi/trustregion
   :alt: Latest PyPI version

This package provides Python routines for solving the trust-region subproblem from nonlinear, nonconvex optimization. For more details on trust-region methods, see the book: A. R. Conn, N. I. M. Gould and Ph. L. Toint (2000), Trust-Region Methods, MPS-SIAM Series on Optimization.

The trust-region subproblem we solve is

.. code-block::

   min_{s in R^n}  g^T s + 0.5 s^T H s, subject to ||s||_2 <= delta (and sl <= s <= su)

**Quick install**

 .. code-block:: bash

    $ sudo apt-get install gfortran
    $ pip install --user numpy
    $ pip install --user trustregion

For more details, see below. Note that NumPy must be installed first, as it is used to compile the Fortran-linked modules.

**Interface** 

The Python package :code:`trustregion` provides one routine, :code:`solve`, with interface:

 .. code-block:: python

    import trustregion
    s               = trustregion.solve(g, H, delta, sl=None, su=None, verbose_output=False)
    s, gnew, crvmin = trustregion.solve(g, H, delta, sl=None, su=None, verbose_output=True)

where the inputs are

* :code:`g`, the gradient of the objective (as a 1D NumPy array)
* :code:`H`, the symmetric Hessian matrix of the objective (as a 2D square NumPy array) - this can be :code:`None` if the model is linear
* :code:`delta`, the trust-region radius (non-negative float)
* :code:`sl`, the lower bounds on the step (as a 1D NumPy array) - this can be :code:`None` if not present, but :code:`sl` and :code:`su` must either be both :code:`None` or both set
* :code:`su`, the upper bounds on the step (as a 1D NumPy array) - this can be :code:`None` if not present, but :code:`sl` and :code:`su` must either be both :code:`None` or both set
* :code:`verbose_output`, a flag indicating which outputs to return.

The outputs are:

* :code:`s`, an approximate minimizer of the subproblem (as a 1D NumPy array)
* :code:`gnew`, the gradient of the objective at the solution :code:`s` (i.e. :code:`gnew = g + H.dot(s)`)
* :code:`crvmin`, a float giving information about the curvature of the problem. If :code:`s` is on the trust-region boundary (given by :code:`delta`), then :code:`crvmin=0`. If :code:`s` is constrained in all directions by the box constraints, then :code:`crvmin=-1`. Otherwise, :code:`crvmin>0` is the smallest curvature seen in the Hessian.

**Example Usage** 

Examples for the use of :code:`trustregion.solve` can be found in the `examples <https://github.com/lindonroberts/trust-region/tree/master/examples>`_ directory on Github.

**Algorithms**

:code:`trustregion` implements three different methods for solving the subproblem, based on the problem class (in Fortran 90, wrapped to Python):

* :code:`trslin.f90` solves the linear objective case (where :code:`H=None` or :code:`H=0`), using Algorithm B.1 from: L. Roberts (2019), `Derivative-Free Algorithms for Nonlinear Optimisation Problems <https://ora.ox.ac.uk/objects/uuid:ec76e895-6eee-491a-88ed-b4ed10fa6003>`_, PhD Thesis, University of Oxford.
* :code:`trsapp.f90` solves the quadratic case without box constraints. It is a minor modification of the routine of the same name in :code:`NEWUOA` [M. J. D. Powell (2004), `The NEWUOA software for unconstrained optimization without derivatives <http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2004_08.pdf>`_, technical report DAMTP 2004/NA05, University of Cambridge].
* :code:`trsbox.f90` solves the quadratic case with box constraints. It is a minor modification of the routine of the same name in :code:`BOBYQA` [M. J. D. Powell (2009), `The BOBYQA algorithm for bound constrained optimization without derivatives <http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf>`_, technical report DAMTP 2009/NA06, University of Cambridge].

In the linear case, an active-set method is used to solve the resulting convex problem. In the quadratic cases, a modification of the Steihaug-Toint/conjugate gradient method is used. For more details, see the relevant references above.

Requirements
------------
:code:`trustregion` requires the following software to be installed:

* Fortran compiler (e.g. gfortran)
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
