import setuptools

# Get package version without "import trustregion" (which requires dependencies to already be installed)
import os
version = {}
with open(os.path.join('trustregion', 'version.py')) as fp:
    exec(fp.read(), version)
__version__ = version['__version__']

try:
    from numpy.distutils.core import setup, Extension
except ModuleNotFoundError:
    print("Please install NumPy before installing trustregion (required to build Fortran subroutines)")
    raise

# Build Fortran module using NumPy (f2Py)
# https://stackoverflow.com/questions/55352409/python-setuptools-compile-fortran-code-and-make-an-entry-points
src_folder = 'trustregion'
sources = []
for f in os.listdir(src_folder):
    if f.endswith('.f90'):
        sources.append(os.path.join(src_folder, f))
lib = Extension(name='trustregion._trs', sources=sources)

# Main package setup
setup(
    name='trustregion',
    version=__version__,
    description='Trust-region subproblem solvers for nonlinear/nonconvex optimization',
    long_description=open('README.rst').read(),
    author='Lindon Roberts',
    author_email='lindon.roberts@anu.edu.au',
    url="https://github.com/lindonroberts/trust-region",
    download_url="https://github.com/lindonroberts/trust-region/archive/v1.2.tar.gz",
    packages=['trustregion'],
    ext_modules=[lib],
    license='GNU GPL',
    keywords = "mathematics nonlinear nonconvex optimization trust-region",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Financial and Insurance Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],
    install_requires = ['numpy'],
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe = True,
    )