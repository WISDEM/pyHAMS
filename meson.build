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
ff = meson.get_compiler('fortran')

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
_ldflags = ['-fopenmp']
if ff.get_id() == 'gcc'
  _fflags = ff.get_supported_arguments('-fdec-math')
  _fflags += ff.get_supported_arguments('-fno-align-commons')
  #_ldflags += ['-llapack']
endif
if ff.get_id() == 'intel'
  _fflags = ff.get_supported_arguments('-mkl')
  #_ldflags = ['-lmkl_rt']
elif ff.get_id() == 'intel-cl'
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _fflags = ff.get_supported_arguments('/mkl')
  #_ldflags += ['-lmkl_rt']
endif
add_project_arguments(_fflags, language: 'fortran')
add_project_link_arguments(_ldflags, language: ['fortran'])

# Note that this works as long as BLAS and LAPACK are detected properly via
# pkg-config. By default we look for OpenBLAS, other libraries can be configured via
# `meson configure -Dblas=blas -Dlapack=lapack` (example to build with Netlib
# BLAS and LAPACK).
# For MKL and for auto-detecting one of multiple libs, we'll need a custom
# dependency in Meson (like is done for scalapack) - see
# https://github.com/mesonbuild/meson/issues/2835
#blas_name = get_option('blas')
#lapack_name = get_option('lapack')
# pkg-config uses a lower-case name while CMake uses a capitalized name, so try
# that too to make the fallback detection with CMake work
#if blas_name == 'openblas'
#  blas_name = ['openblas', 'OpenBLAS']
#endif
#if lapack_name == 'openblas'
#  lapack_name = ['openblas', 'OpenBLAS']
#endif
#blas = dependency(blas_name)
lapack = dependency('lapack')

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
