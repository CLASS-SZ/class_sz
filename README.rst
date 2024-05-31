==============================================
CLASS_SZ
==============================================
 Cosmic Linear Anisotropy Solving System with Fast and Accurate CMB, LSS and Halo Model Observables Computations

To install/run CLASS_SZ on the fly, check the colab notebook: 

`class_sz colab notebook <https://colab.research.google.com/drive/1AULgG4ZLLG1YXRI86L54-hpjWyl1X-8c?usp=sharing>`_.

Further informaton on dowloading and installing the code are given hereafter, for Mac, M1/2 Mac, and Linux.

The tutorial notebooks can be found at:

https://github.com/CLASS-SZ/notebooks

These notebooks along with the paper (`Bolliet et al 2023 <https://arxiv.org/abs/2310.18482>`_) constitute the current documentation.

CLASS_SZ is as fast as it gets, with full parallelization, implementation of high-accuracy cosmopower emulators (see below for some instructions) and usage of Fast Fourier Transforms (including FFTLog).

Since it is based on Lesgourgues's class code, the halo model and LSS calculations (essentially based on distances and
matter clustering) are always consistent with the cosmological model computed by class.

CLASS_SZ is an extension of Julien Lesgourgues's `CLASS <https://github.com/lesgourg/class_public>`_ code.

CLASS_SZ is initially based on Eiichiro Komatsu’s fortran code `SZFAST <http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/>`_.

CLASS_SZ functionalities are located in the files:

- **source/class_sz.c** for the main CLASS_SZ functions, 

- **tools/class_sz_tools.c** for useful routines,

- **source/class_sz_clustercounts.c** for tSZ cluster counts. Since March 2024, CLASS_SZ cluster counts calculations are superseded by `cosmocnc <https://github.com/inigozubeldia/cosmocnc>`_ (`Zubeldia & Bolliet 2024 <https://arxiv.org/abs/2403.09589>`_). 

CLASS_SZ's outputs are regularly cross-checked with other CMBxLSS codes, such as:

- `cosmocnc <https://github.com/inigozubeldia/cosmocnc>`_,

- `hmvec <https://github.com/simonsobs/hmvec/tree/master/hmvec>`_,

- `ccl <https://github.com/LSSTDESC/CCL>`_,

- `HaloGen <https://github.com/EmmanuelSchaan/HaloGen/tree/master>`_,

- `yxg <https://github.com/nikfilippas/yxg>`_,

- `halomodel_cib_tsz_cibxtsz <https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz>`_.





Downloading the code
--------------------

Clone or download from https://github.com/CLASS-SZ/class_sz

Note: the significant size of the repository is due to the size of the original **class** repository.


Using the code
--------------

The **class_sz** code is public.


If you use it, please cite:

`CLASS_SZ: I Overview (Boris Bolliet, Aleksandra Kusiak, Fiona McCarthy, et al, 2024) <https://arxiv.org/abs/2310.18482>`_.

`Projected-field kinetic Sunyaev-Zel'dovich Cross-correlations: halo model and forecasts (Boris Bolliet, J. Colin Hill, Simone Ferraro, Aleksandra Kusiak, Alex Krolewski, 2023) <https://iopscience.iop.org/article/10.1088/1475-7516/2023/03/039>`_.

If you use the emulators (fast method of class_sz, see below), please cite:

`High-accuracy emulators for observables in LCDM, Neff+LCDM, Mnu+LCDM and wCDM cosmologies (Bolliet, Spurio Mancini, Hill, Madhavacheril, Jense, Calabrese, Dunkley, 2023) <https://inspirehep.net/literature/2638458>`_.

`COSMOPOWER: emulating cosmological power spectra for accelerated Bayesian inference from next-generation surveys (Spurio Mancini, Piras, Alsing, Joachimi, Hobson, 2021) <https://arxiv.org/abs/2106.03846>`_.


If you use thermal SZ power spectrum and cluster counts calculations, please also consider citing:

`Including massive neutrinos in thermal Sunyaev Zeldovich power spectrum and cluster counts analyses (Bolliet, Brinckmann, Chluba, Lesgourgues, 2020) <https://arxiv.org/abs/1906.10359>`_.

`Dark Energy from the Thermal Sunyaev Zeldovich Power Spectrum (Bolliet, Comis, Komatsu, Macias-Perez, 2017)
<https://arxiv.org/abs/1712.00788>`_.

`The Sunyaev-Zel'dovich angular power spectrum as a probe of cosmological parameters (Komatsu and Seljak, 2002)
<https://arxiv.org/abs/astro-ph/0205468>`_.

If you use the code, please also cite the original class papers (since class_sz is an extension of class), e.g.,:

`CLASS I: Overview (Lesgourgues, 2011) <https://arxiv.org/abs/1104.2932>`_.

`CLASS II: Approximation schemes (Blas, Lesgourgues, Tram, 2011)
<http://arxiv.org/abs/1104.2933>`_.

As well as other references listed there: http://class-code.net



Compiling CLASS_SZ and Getting Started
======================================

A Colab notebook shows you how to quick start. You can even run your calculations there (although this may not be as fast as on your computer since Colab, as of 2023, runs on two cores):

`class_sz Colab notebook <https://colab.research.google.com/drive/1AULgG4ZLLG1YXRI86L54-hpjWyl1X-8c?usp=sharing>`_.

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

    echo -e "\n# Set library paths for class_sz\nexport LIBRARY_PATH=/Users/boris/opt/miniconda3/lib:path/to/gsl/:path/to/fftw/:$LIBRARY_PATH\nexport C_INCLUDE_PATH=/Users/boris/opt/miniconda3/include/:path/to/gsl/:path/to/fftw/:$C_INCLUDE_PATH\nexport DYLD_LIBRARY_PATH=\"/Users/boris/opt/miniconda3/lib:\$DYLD_LIBRARY_PATH\" # (Mac M1 users only)" >> ~/.bashrc

To apply the changes immediately:

.. code-block:: bash

    source ~/.bashrc

For `~/.bash_profile` (common for macOS):

.. code-block:: bash

    echo -e "\n# Set library paths for class_sz\nexport LIBRARY_PATH=/Users/boris/opt/miniconda3/lib:path/to/gsl/:path/to/fftw/:$LIBRARY_PATH\nexport C_INCLUDE_PATH=/Users/boris/opt/miniconda3/include/:path/to/gsl/:path/to/fftw/:$C_INCLUDE_PATH\nexport DYLD_LIBRARY_PATH=\"/Users/boris/opt/miniconda3/lib:\$DYLD_LIBRARY_PATH\" # (Mac M1 users only)" >> ~/.bash_profile

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


Computing CMB, LSS and halo model quantities via the Python wrapper classy_sz
-----------------------------------------------------------------------------

Class_sz is now very fast ! In part it's because it can run with emulators. This is available via the python wrapper (if requested).

Once class_sz is installed. You can use the python wrapper classy_sz just as you use classy with the normal class code.
You can compute everything classy computes, as well as all the additional CMB, LSS and Halo Model quantities implemented in class_sz.

First, make sure that you have compiled the python wrapper with:

$ make clean

$ make -j

(Note that the second command must be 'make -j', and not 'make -j class_sz' for the python wrappper to be compiled.)

Have a look at the notebooks https://github.com/CLASS-SZ/notebooks. They all use the python wrapper.


Python Wrapper (Tensorflow and Cosmopower Dependency)
-----------------------------------------------------

Since recently we have implemented emulators in classy_sz, now it has an extra-dependency to tensorflow through cosmopower.

1. Install tensoflow first (see below for Mac M1 specific issues).
2. Then install cosmopower (https://alessiospuriomancini.github.io/cosmopower/installation/). Note that the needed tensorflow version may not be the lattest, see the requirements (https://github.com/alessiospuriomancini/cosmopower/blob/main/requirements.txt). 
3. Clone the https://github.com/cosmopower-organization/notebooks repo.
4. Open notebooks/get_quantities_cosmopower.ipynb notebook and follow the instructions there to get the cosmopower emulators.
5. $ make -j
6. Check you can import classy_sz in your python/jupyter notebook, e.g.,:
  $ python

  $ from classy_sz import Class
or try to run any of the notebooks.

7. To run the emulator-based computations, simply change
  M.compute()

to

  M.compute_class_szfast()

8. There are many examples in the notebooks how to use class_szfast. See https://github.com/CLASS-SZ/notebooks.




Some tips to run on computer clusters
---------------------------------------

Module load, module show to get gsl and fftw.
At NERC/Cori, the code works with gsl/2.7. (There seems to be a problematic behavior during job submission with gsl/2.5.)

Mpi4py needs to be correctly installed. Follow:
https://cobaya.readthedocs.io/en/latest/installation.html#mpi-parallelization-optional-but-encouraged
You may need to activate an environment to run the install comment.
To make sure you use the same openmpi compiler, example:
env MPICC=/global/common/software/m3169/cori/openmpi/4.1.2/intel/bin/mpicc python -m pip install mpi4py



GSL library
------------------------------


New version of class_sz requires gsl (for the integration routines).
One may need to edit the **Makefile** adding the include path for gsl libraries, e.g.,:


    INCLUDES = -I../include -I/usr/local/include/ **-I/path_to_gsl/gsl-2.6/include/**

    class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS) $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -lm **-L/path_to_gsl/gsl-2.6/lib/ -lgsl -lgslcblas** -lfftw3

For the python wrapper, one also may need to add the absolute path to gsl libraries, e.g.,:

in **class_sz/python/setup.py**:

    classy_ext = Extension("classy", [os.path.join(classy_folder, "classy.pyx")], include_dirs=[nm.get_include(), include_folder, '**/path/to/gsl-2.6/include**'], libraries=liblist,library_dirs=[root_folder, GCCPATH],extra_link_args=['-lgomp','**-L/path_to_gsl/gsl-2.6/lib/**','**-lgsl**','**-lgslcblas**',-lfftw3])



When running, the gsl library also need to be included in the environment variables, i.e., one may
need to do:

    $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path_to_gsl/gsl-2.6/lib

    $ export LD_LIBRARY_PATH

Note that these prescriptions are system dependent: you may not need them if your path and environment variables are such that gsl and its libraries are well linked.
If you are tired of having to execute these lines each time you run codes in a fresh terminal, just paste them in your bash profile file (the one that ends with .sh).

FFTLog library
------------------------------

class_sz now requires FFTW3 library.

If the code complains about the library not being found, just make sure you followed the same installation instruction as you did for gsl.
Namely, edit the the Makefile with the path to the include files (the ones that end with '.h') -I/path_to_fftw3/fftw3/include/, the path to the library files (the ones that end with .so,.a, .dylib, and so on) -L/path_to_fftw3/fftw3/lib/. The setup.py file may also need to be amended accordingly.
And also make sure you do:

    $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path_to_fftw3/fftw3/lib

    $ export LD_LIBRARY_PATH

if the previous modifs were not enough.

MacOS problem with OpenMP
------------------------------

To run the code in parallel, you may run into a problem on a mac. The solution is provided here:

https://github.com/lesgourg/class_public/issues/208

Essentially, you need to edit a line in python/setup.py such as the code knows about the mpi libraries to be used with your compiler (gcc-11 in the example below).
In our case the modif looks like this:

  extra_link_args=['-lgomp','-lgsl','-lgslcblas','**-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/11/**']

Mac OS with M1 chip
----------------------

We advise installing fftw, gsl, openmp with anaconda, i.e., conda forge etc..

LD_LIBRARY_PATH becomes DYLD_LIBRARY_PATH, hence, export with, e.g.,:

export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/anaconda3/lib


In Makefile:

CC = clang
PYTHON ?= /set/path/to/anaconda3/python
OPTFLAG = -O4 -ffast-math # dont use: -arch x86_64
OMPFLAG   = -Xclang -fopenmp
LDFLAG += -lomp
INCLUDES =  -I../include -I/usr/local/include/ -I/path/to/anaconda3/include/
$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -L/usr/local/lib -L/path/to/anaconda3/lib/ -lgsl -lgslcblas -lfftw3 -lm

In setup.py:

extra_link_args=['-lomp','-lgsl','-lfftw3','-lgslcblas'])

Tensorflow on mac M1
----------------------

To install the new version of class_sz, you will need tensorflow (needed for the cosmopower emulators). On M1/M2 make sure, you have the arch64 version of conda (if not, you need to remove your entire conda and install the arch64 version for Apple sillicon).

This video might be helpful https://www.youtube.com/watch?v=BEUU-icPg78
Then you can follow standard Tensorflow installation recipe for M1, e.g., https://caffeinedev.medium.com/how-to-install-tensorflow-on-m1-mac-8e9b91d93706 or https://developer.apple.com/forums/thread/697846 .

Compiler - GCC version
------------------------------

The gcc copiler can be changed easily to any gcc version that is available to you.
There are two modifications:

1) Line 20 of Makefile: CC = gcc-XX (where XX=11 in our case.)

2) Line 12 of python/setup.py: replace 'gcc-11' with, e.g., 'gcc-XX'.  


Pre M1 Mac  
----------------------  
See Makefile_preM1mac for an example makefile for older Macs (without the M1 chip). Some key points include adding paths involving libomp to LDFLAG and INCLUDES.
In python/setup.py, you may also want to modify the extra_link_args list to contain '-lomp' instead of '-lgomp' and add the libomp library path as well to that list. 
For example, extra_link_args=['-lomp', '-lgsl','-lfftw3','-lgslcblas', '-L/usr/local/opt/libomp/lib/'].  





Support
-------

To get support on the class_sz module, feel free to open an issue on the GitHub page, we will try to answer as soon as possible.
