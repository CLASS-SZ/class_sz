==============================================
CLASS_SZ
==============================================
Cosmic Linear Anisotropy Solving System with Machine Learning Accelerated and Accurate CMB, LSS, and Halo Model Observables Computations

To install/run CLASS_SZ on the fly, check the Colab notebook:

`Colab Notebook <https://colab.research.google.com/drive/1AULgG4ZLLG1YXRI86L54-hpjWyl1X-8c?usp=sharing>`_.

Further information on downloading and installing the code is given hereafter for Mac, M1/2 Mac, and Linux.

The tutorial notebooks can be found at:

`CLASS_SZ Notebooks <https://github.com/CLASS-SZ/notebooks>`_

These notebooks along with the paper (`Bolliet et al 2023 <https://arxiv.org/abs/2310.18482>`_) constitute the current documentation.

CLASS_SZ is as fast as it gets, with full parallelization, implementation of high-accuracy neural network emulators, and Fast Fourier Transforms.

Since it is based on Lesgourgues's class code, the halo model and LSS calculations (essentially based on distances and matter clustering) are always consistent with the cosmological model computed by CLASS.

CLASS_SZ is an extension of Julien Lesgourgues's `CLASS <https://github.com/lesgourg/class_public>`_ code.

CLASS_SZ is initially based on Eiichiro Komatsu’s Fortran code `SZFAST <http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/>`_.

CLASS_SZ functionalities are located in the files:

- **source/class_sz.c** for the main CLASS_SZ functions,
- **tools/class_sz_tools.c** for useful routines,
- **source/class_sz_clustercounts.c** for tSZ cluster counts. Since March 2024, CLASS_SZ cluster counts calculations are superseded by `cosmocnc <https://github.com/inigozubeldia/cosmocnc>`_ (`Zubeldia & Bolliet 2024 <https://arxiv.org/abs/2403.09589>`_).

CLASS_SZ's outputs are regularly cross-checked with other CMBxLSS codes, such as:

- `cosmocnc <https://github.com/inigozubeldia/cosmocnc>`_
- `hmvec <https://github.com/simonsobs/hmvec/tree/master/hmvec>`_
- `ccl <https://github.com/LSSTDESC/CCL>`_
- `HaloGen <https://github.com/EmmanuelSchaan/HaloGen/tree/master>`_
- `yxg <https://github.com/nikfilippas/yxg>`_
- `halomodel_cib_tsz_cibxtsz <https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz>`_

Downloading the Code
--------------------

Clone or download from `GitHub <https://github.com/CLASS-SZ/class_sz>`_

Note: the significant size of the repository is due to the size of the original **class** repository.

Using the Code
--------------

The **class_sz** code is public.

If you use it, please cite:

- `CLASS_SZ: I Overview (Bolliet et al. 2024) <https://arxiv.org/abs/2310.18482>`_
- `Projected-field kinetic Sunyaev-Zel'dovich Cross-correlations: halo model and forecasts (Bolliet et al. 2023) <https://iopscience.iop.org/article/10.1088/1475-7516/2023/03/039>`_

If you use the emulators (fast method of class_sz, see below), please cite:

- `High-accuracy emulators for observables in LCDM, Neff+LCDM, Mnu+LCDM and wCDM cosmologies (Bolliet et al. 2023) <https://inspirehep.net/literature/2638458>`_
- `COSMOPOWER: emulating cosmological power spectra for accelerated Bayesian inference from next-generation surveys (Spurio Mancini et al. 2021) <https://arxiv.org/abs/2106.03846>`_

If you use thermal SZ power spectrum and cluster counts calculations, please also consider citing:

- `Including massive neutrinos in thermal Sunyaev Zeldovich power spectrum and cluster counts analyses (Bolliet et al. 2020) <https://arxiv.org/abs/1906.10359>`_
- `Dark Energy from the Thermal Sunyaev Zeldovich Power Spectrum (Bolliet et al. 2017) <https://arxiv.org/abs/1712.00788>`_
- `The Sunyaev-Zel'dovich angular power spectrum as a probe of cosmological parameters (Komatsu and Seljak, 2002) <https://arxiv.org/abs/astro-ph/0205468>`_

If you use the code, please also cite the original CLASS papers (since class_sz is an extension of CLASS), e.g.,:

- `CLASS I: Overview (Lesgourgues, 2011) <https://arxiv.org/abs/1104.2932>`_
- `CLASS II: Approximation schemes (Blas, Lesgourgues, Tram, 2011) <http://arxiv.org/abs/1104.2933>`_

As well as other references listed here: http://class-code.net

Compiling CLASS_SZ and Getting Started
======================================

A Colab notebook shows you how to quick start. You can even run your calculations there (although this may not be as fast as on your computer since Colab, as of 2023, runs on two cores):

`Colab Notebook <https://colab.research.google.com/drive/1AULgG4ZLLG1YXRI86L54-hpjWyl1X-8c?usp=sharing>`_.

Move to the code repository:

.. code-block:: bash

    $ cd class_sz

Pick the correct Makefile:

For Mac M1:

.. code-block:: bash

    $ cp Makefile_m1 Makefile

For Linux:

.. code-block:: bash

    $ cp Makefile_linux Makefile

Clean up and compile:

.. code-block:: bash

    $ make clean
    $ make -j

(You may need to use `$ sudo make`.)

The previous commands compile both the executable and the Python wrapper. The `-j` flag speeds up the compilation process by using multiple cores.

For Mac users, class_sz also works on the Mac M1 chips. M2 chips have not been tested yet.

Library Path Configuration
--------------------------

It is often the case that some libraries are not found. In general, setting the following paths appropriately should solve your issues:

.. code-block:: bash

    export LIBRARY_PATH=/Users/boris/opt/miniconda3/lib:path/to/gsl/:path/to/fftw/:$LIBRARY_PATH
    export C_INCLUDE_PATH=/Users/boris/opt/miniconda3/include/:path/to/gsl/:path/to/fftw/:$C_INCLUDE_PATH
    export DYLD_LIBRARY_PATH="/Users/boris/opt/miniconda3/lib:$DYLD_LIBRARY_PATH" # (Mac M1 users only)

To ensure these paths are set every time you open a terminal, you can add these lines to your `~/.bashrc` or `~/.bash_profile` file automatically using the `echo` command.

For `~/.bashrc` (common for most Linux systems):

.. code-block:: bash

    echo -e "\n# Set library paths for class_sz\nexport LIBRARY_PATH=/Users/boris/opt/miniconda3/lib:path/to/gsl/:path/to/fftw/:\$LIBRARY_PATH\nexport C_INCLUDE_PATH=/Users/boris/opt/miniconda3/include/:path/to/gsl/:path/to/fftw/:\$C_INCLUDE_PATH\nexport DYLD_LIBRARY_PATH=\"/Users/boris/opt/miniconda3/lib:\$DYLD_LIBRARY_PATH\" # (Mac M1 users only)" >> ~/.bashrc

To apply the changes immediately:

.. code-block:: bash

    source ~/.bashrc

For `~/.bash_profile` (common for macOS):

.. code-block:: bash

    echo -e "\n# Set library paths for class_sz\nexport LIBRARY_PATH=/Users/boris/opt/miniconda3/lib:path/to/gsl/:path/to/fftw/:\$LIBRARY_PATH\nexport C_INCLUDE_PATH=/Users/boris/opt/miniconda3/include/:path/to/gsl/:path/to/fftw/:\$C_INCLUDE_PATH\nexport DYLD_LIBRARY_PATH=\"/Users/boris/opt/miniconda3/lib:\$DYLD_LIBRARY_PATH\" # (Mac M1 users only)" >> ~/.bash_profile

To apply the changes immediately:

.. code-block:: bash

    source ~/.bash_profile

Running the Code
----------------

Run the code with most of the power spectra output:

.. code-block:: bash

    $ ./class_sz class_sz_test.ini

The `.ini` files are the parameter files.

If you want to run CLASS and not do the class_sz part, you can! For example:

.. code-block:: bash

    $ ./class_sz explanatory.ini

This will just run the standard CLASS code and its calculations. All depends on what output you request: if you request a class_sz observable or not.

Computing CMB, LSS and Halo Model Quantities via the Python Wrapper classy_sz
-----------------------------------------------------------------------------

CLASS_SZ is now very fast! In part it's because it can run with emulators. This is available via the Python wrapper (if requested).

Once CLASS_SZ is installed, you can use the Python wrapper classy_sz just as you use classy with the normal CLASS code. You can compute everything classy computes, as well as all the additional CMB, LSS, and Halo Model quantities implemented in CLASS_SZ.

First, make sure that you have compiled the Python wrapper with:

.. code-block:: bash

    $ make clean
    $ make -j

(Note that the second command must be `$ make -j`, and not `$ make -j class_sz` for the Python wrapper to be compiled.)

Have a look at the notebooks `CLASS_SZ Notebooks <https://github.com/CLASS-SZ/notebooks>`_. They all use the Python wrapper.

Accelerated Computations
------------------------

To run the emulator-based computations, simply change:

.. code-block:: python

    M.compute()

to:

.. code-block:: python

    M.compute_class_szfast()

There are many examples in the notebooks on how to use class_szfast. See `CLASS_SZ Notebooks <https://github.com/CLASS-SZ/notebooks>`_.

Some Tips to Run on Computer Clusters
-------------------------------------

Use module load, module show to get GSL and FFTW.
At NERSC/Cori/Perlmutter, the code works with gsl/2.7. (There seems to be a problematic behavior during job submission with gsl/2.5.)

For Monte Carlo analyses, we also recall that Mpi4py needs to be correctly installed. Follow:
https://cobaya.readthedocs.io/en/latest/installation.html#mpi-parallelization-optional-but-encouraged

TensorFlow on Mac M1
--------------------

To install the new version of CLASS_SZ, you will need TensorFlow (needed for the Cosmopower emulators). On M1/M2, make sure you have the arm64 version of conda (if not, you need to remove your entire conda and install the arm64 version for Apple Silicon).

This video might be helpful: `Installing TensorFlow on M1 Mac <https://www.youtube.com/watch?v=BEUU-icPg78>`_.
Then you can follow the standard TensorFlow installation recipe for M1, e.g., `Medium Article <https://caffeinedev.medium.com/how-to-install-tensorflow-on-m1-mac-8e9b91d93706>`_ or the `Apple Developer Forums <https://developer.apple.com/forums/thread/697846>`_.

The following two lines should fix most issues:

.. code-block:: bash

    pip install tensorflow-metal-1.1.0
    conda install -c apple tensorflow-deps

Pre M1 Mac
----------

See Makefile_preM1mac for an example makefile for older Macs (without the M1 chip). Some key points include adding paths involving libomp to LDFLAG and INCLUDES.
In python/setup.py, you may also want to modify the extra_link_args list to contain '-lomp' instead of '-lgomp' and add the libomp library path as well to that list. 
For example, extra_link_args=['-lomp', '-lgsl','-lfftw3','-lgslcblas', '-L/usr/local/opt/libomp/lib/'].

This makefile is not maintained anymore but we keep it for reference. If you need to run class_sz on a pre-M1 Mac and have serious issues, please contact us.

Support
-------

To get support on the class_sz module, feel free to open an issue on the GitHub page, we will try to answer as soon as possible.
