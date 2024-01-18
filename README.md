# pyHAMS

This is a python wrapper for [HAMS](https://github.com/YingyiLiu/HAMS), a boundary element method for hydrodynamic analysis of submerged structures. 

There is cylinder test case that demonstrates usage and outputs in ``test/test_cylinder.py``.

## Prerequisites

pyHAMS requires a Fortran compiler.  The python wrapper currently supports GNU and Intel compilers.  HAMS can be built with Flang, but that is not yet recognized by pyHAMS.

## Install (as a library)

To install pyHAMS as a library that can be used by WEIS or RAFT in the backend, you can use either conda or PyPi package management:

    $ pip install pyHAMS
	
or
	
    $ conda install pyHAMS


## Install (from source)

If you would like to build the project locally from source for easier access to the underlying methods and tests, do:

    $ git clone https://github.com/WISDEM/pyHAMS.git
    $ cd pyHAMS
    $ pip install .

If developer/editable mode, do the same `git clone` step, but on install do:

    $ pip install --no-build-isolation -e .

The `--no-build-isolation` option is important per [Meson guidelines](https://meson-python.readthedocs.io/en/latest/how-to-guides/editable-installs.html).  Note that this package uses `mesonpy` for an installation backend, so there is no `setup.py` file.


## Unit Tests

    $ pytest test

For software issues please use <https://github.com/WISDEM/pyHAMS/issues>.  For theory related questions and comments please use the [HAMS issues page](https://github.com/YingyiLiu/HAMS/issues).


## License

pyHAMS uses the Apache License v2.0

