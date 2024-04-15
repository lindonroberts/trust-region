incdir_numpy = run_command(py3, 
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'], 
    check: true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

py3.extension_module('_trs',
    ['trsapp.f90', 'trsbox.f90', 'trslin.f90', '_trsmodule.c'],
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    dependencies: py3_dep,
    install: true,
    subdir: 'trustregion'
)

python_sources = ['__init__.py', 'interface.py']

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'trustregion'
)
