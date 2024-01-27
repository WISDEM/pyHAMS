# pyHAMS

This is a python wrapper for [HAMS](https://github.com/YingyiLiu/HAMS), a boundary element method for hydrodynamic analysis of submerged structures. 

There is cylinder test case that demonstrates usage and outputs in ``test/test_cylinder.py``.

## Prerequisites

pyHAMS requires a Fortran compiler and OpenBLAS / MKL / LAPACK.  We strongly recommend the [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3) distribution of `conda` to satisfy package dependencies, as the traditional Anaconda distribution can be slow and struggle with mapping out dependencies.

## Install (as a library)

To install pyHAMS as a library that can be used by WEIS or RAFT in the backend, conda is your best option:
	
    $ conda install pyHAMS


## Install (from source)

If you would like to build the project locally from source for easier access to the underlying methods and tests, we still recommend using `conda` to satisfy dependencies.

    $ git clone https://github.com/WISDEM/pyHAMS.git
    $ cd pyHAMS
    $ conda env create --name pyhams-env -f environment.yml
    $ conda activate pyhams-env
    $ conda install -y compilers                       # (Mac/Linux without other compilers)
    $ conda install -y m2w64-toolchain libpython       # (Windows only)
    $ pip install .

If developer/editable mode, replace the final step with:

    $ pip install --no-build-isolation -e .

The `--no-build-isolation` option is important per [Meson guidelines](https://meson-python.readthedocs.io/en/latest/how-to-guides/editable-installs.html).  Note that this package uses `mesonpy` for an installation backend, so there is no `setup.py` file.


## Unit Tests

    $ pytest test

For software issues please use <https://github.com/WISDEM/pyHAMS/issues>.  For theory related questions and comments please use the [HAMS issues page](https://github.com/YingyiLiu/HAMS/issues).


## License

pyHAMS uses the Apache License v2.0

