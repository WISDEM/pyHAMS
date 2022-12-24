# Much of this is from SciPy

project(
  'pyhams',
  'c',
# unnecessary metadata commented out until Meson supports PEP517 and installation with pip
  license: 'Apache',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c11',
  ],
)

cc = meson.get_compiler('c')
add_languages('fortran', native: false)
fc = meson.get_compiler('fortran')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

is_windows = host_machine.system() == 'windows'

# Adding at project level causes many spurious -lgfortran flags.
_fflags = []
if fc.get_id() == 'gcc'
  _fflags = fc.get_supported_arguments('-fdec-math')
  _fflags += fc.get_supported_arguments('-fno-align-commons')
endif
if fc.get_id() == 'intel'
  _fflags = fc.get_supported_arguments('-mkl')
elif fc.get_id() == 'intel-cl'
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _fflags = fc.get_supported_arguments('/mkl')
endif
add_project_arguments(_fflags, language: 'fortran')
omp = dependency('openmp', required: false)
if not omp.found()
  if fc.get_id() == 'intel'
    omp = declare_dependency(
      link_args: '-qopenmp',
    )
  else
    omp = declare_dependency(
      link_args: '-fopenmp',
    )
  endif
endif

if fc.get_id() == 'gcc'
  lapack = dependency('openblas', required: false)
  if not lapack.found()
    lapack = fc.find_library('openblas', required: false)
  endif
  if not lapack.found()
    lapack = dependency('lapack', required: false)
  endif
  if not lapack.found()
    lapack = fc.find_library('lapack', required: false)
  endif
  if not lapack.found()
    lapack = declare_dependency(
      link_args: '-llapack',
    )
  endif
endif
if fc.get_id() == 'intel'
  lapack = dependency('mkl', required: false)
  if not lapack.found()
    lapack = fc.find_library('mkl_rt')
  endif
  if not lapack.found()
    lapack = declare_dependency(
      link_args: '-lmkl_rt',
    )
  endif
endif

# https://mesonbuild.com/Python-module.html
# Here we differentiate from the python used by meson, py3_command, and that python target, py3_target. This is useful
# when cross compiling like on conda-forge
py_mod = import('python')
py3 = py_mod.find_installation('python')
py3_dep = py3.dependency()

message(py3.path())
message(py3.get_install_dir())

if not is_windows
  subdir('pyhams')
endif