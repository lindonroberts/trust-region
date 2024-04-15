# When updating versions, 
# 1. Increment version number both here and in trustregion/__init__.py 
# 2. In trustregion directory, run
#        $ python -m numpy.f2py *.f90 -m _trs
#    to enable the build to work propertly

project('trustregion', 'c', version: '1.2.1', default_options: ['c_std=c99', 'fortran_std=legacy'])

add_languages('fortran', native: false)

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()

subdir('trustregion')

