# NumPy include directory - needed in all submodules
incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(py3,
                             [
                               '-c',
                               'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
                             ],
                             check: true
                            ).stdout().strip()
endif
# this creates a raw string which is useful for Windows use of '\' for paths
incdir_numpy = '''@0@'''.format(incdir_numpy)
inc_np = include_directories(incdir_numpy)
numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']
np_dep = declare_dependency(include_directories: inc_np, compile_args: numpy_nodepr_api)

#incdir_f2py = run_command(py3,
#    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
#    check : true
#).stdout().strip()

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
fortranobject_c = incdir_f2py / 'fortranobject.c'
inc_f2py = include_directories(incdir_f2py)

# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: numpy_nodepr_api,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

hams_connect = custom_target('_hamsmodule.c',
                             input : ['src/_hams.pyf'],
                             output : ['_hamsmodule.c', '_hams-f2pywrappers2.f90'],
                             command: [py3, '-m', 'numpy.f2py',
			               '@INPUT@', '--lower', '--build-dir', 'pyhams']
                            )

hams_sources = ['src/HAMS_LAPACK.f90',
                'src/WavDynMods.f90',
                'src/PatclVelct.f90',
                'src/BodyIntgr.f90',
                'src/BodyIntgr_irr.f90',
                'src/AssbMatx.f90',    
                'src/AssbMatx_irr.f90',
                'src/SingularIntgr.f90',
                'src/InfGreen_Appr.f90',
                'src/FinGrnExtSubs.f90',
                'src/FinGreen3D.f90',
                'src/CalGreenFunc.f90',
                'src/HydroStatic.f90',
                'src/ImplementSubs.f90',
                'src/InputFiles.f90',
                'src/NormalProcess.f90',
                'src/ReadPanelMesh.f90',
                'src/PotentWavForce.f90',
                'src/PressureElevation.f90',
                'src/PrintOutput.f90',
                'src/SolveMotion.f90',
                'src/WavDynSubs.f90',
                'src/HAMS_Prog.f90',
               ]

py3.extension_module('_hams',
                     hams_sources,
                     hams_connect,
                     dependencies : [fortranobject_dep, omp],
                     subdir: 'pyhams',
		     link_language: 'fortran',
                     install : true)


python_sources = [
  '__init__.py',
  'pyhams.py',
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'pyhams'
)
