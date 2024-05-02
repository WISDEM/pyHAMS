# Much of this is from SciPy

project(
  'pyhams',
  'c',
# unnecessary metadata commented out until Meson supports PEP517 and installation with pip
  license: 'Apache',
  #meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c11',
  ],
)

cc = meson.get_compiler('c')
add_languages('fortran', native: false)
fc = meson.get_compiler('fortran')
message(fc.get_id())
is_windows = host_machine.system() == 'windows'

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
if fc.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

_fflags = []
if fc.get_id() in ['intel', 'intel-llvm']
  _fflags = fc.get_supported_arguments('-qmkl')
  _fflags += fc.get_supported_arguments('-fp-model=strict')
  minus0_arg = ['-assume', 'minus0']
  if fc.has_multi_arguments(minus0_arg)
    _fflags += minus0_arg
  endif
elif fc.get_id() in ['intel-cl', 'intel-llvm-cl']
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _fflags = fc.get_supported_arguments('/qmkl')
  _fflags += fc.get_supported_arguments(
    '/fp:strict', '/MD', '/names:lowercase', '/assume:underscore',
    '/assume:minus0'
    )  
else
  _fflags = fc.get_supported_arguments('-fdec-math')
  _fflags += fc.get_supported_arguments('-fno-align-commons')
endif
add_project_arguments(_fflags, language: 'fortran')

omp = dependency('openmp')
#if fc.get_id() in ['intel', 'intel-llvm']
#  omp = declare_dependency(link_args: '-qopenmp')
#elif fc.get_id() in ['intel-cl', 'intel-llvm-cl']
#  omp = declare_dependency(link_args: '/qopenmp')
#else
#  omp = declare_dependency(link_args: '-fopenmp')
#endif

# Need lapack via OpenBlas preferably
lapack = dependency('openblas', ['mkl', 'mkl_rt', 'mkl-sdl'], 'lapack', 'blas', 'accelerate', required: true)

# https://mesonbuild.com/Python-module.html
# Here we differentiate from the python used by meson, py3_command, and that python target, py3_target. This is useful
# when cross compiling like on conda-forge
py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

message(py3.path())
message(py3.get_install_dir())

subdir('pyhams')
