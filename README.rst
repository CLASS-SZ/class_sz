==============================================
CLASS_SZ
==============================================
 Cosmic Linear Anisotropy Solving System

 with fast and accurate halo model computations



In addition to SZ power spectrum, class_sz can compute cross and auto power spectra for other tracers
in the halo model (kSZ, galaxy, galaxy-lensing, ISW, CMB lensing and CIB).

It has several mass functions implemented, with several possible halo mass definitions and concentration-mass
relations. For galaxy clustering and lensing, class_sz has an implementation of HOD based on the one used by
the DES collaboration.

The code is close to be as fast as it can get, with full parallelization.

Since it is based on Lesgourgues's class code, the halo model (essentially based on distances and
matter clustering) is always consistent with the cosmological model.



**Take a look at the notebook to see what class_sz can do:**

https://github.com/borisbolliet/class_sz/blob/master/notebooks/class_sz_plots_and_tutorial.ipynb

The code is currently in development, don't hesitate to reach out if you would like to use the code and need assistance.

CLASS_SZ is an extension of Julien Lesgourgues's CLASS code.

For download and information on CLASS, see http://class-code.net and https://github.com/lesgourg/class_public

CLASS_SZ is initially based on Eiichiro Komatsu’s fortran code SZFAST.

(See http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/)

CLASS_SZ modules are located in the files **source/class_sz.c** and **source/class_sz_clustercounts.c**.


CLASS_SZ's outputs are regularly cross-checked with other halo model codes, such as:

- `hmvec <https://github.com/simonsobs/hmvec/tree/master/hmvec>`_,

- `ccl <https://github.com/LSSTDESC/CCL>`_,

- `HaloGen <https://github.com/EmmanuelSchaan/HaloGen/tree/master>`_,

- `yxg <https://github.com/nikfilippas/yxg>`_.



Downloading the code
--------------

Clone or download from https://github.com/borisbolliet/class_sz

Note: the significant size of the repository is due to the size of the original **class** repository.


Using the code
--------------

The **class_sz** code is public.

Some References.

The first papers using class_sz were:

`Including massive neutrinos in thermal Sunyaev Zeldovich power spectrum and cluster counts analyses (Bolliet, Brinckmann, Chluba, Lesgourgues, 2020) <https://arxiv.org/abs/1906.10359>`_.

`Dark Energy from the Thermal Sunyaev Zeldovich Power Spectrum (Bolliet, Comis, Komatsu, Macias-Perez, 2017)
<https://arxiv.org/abs/1712.00788>`_.

If you use the code, please also cite the original class papers (since class_sz is an extension of class), e.g.,:

`CLASS II: Approximation schemes (Blas, Lesgourgues, Tram, 2011)
<http://arxiv.org/abs/1104.2933>`_.

As well as the original tSZ power spectrum halo-model paper:

`The Sunyaev-Zel'dovich angular power spectrum as a probe of cosmological parameters (Komatsu and Seljak, 2002)
<https://arxiv.org/abs/astro-ph/0205468>`_.


Compiling CLASS_SZ and getting started
--------------------------------------

Move to the code repository

    $ cd class_sz

Clean up and Compile

    $ make clean

    $ make

(You may need to do a ‘$ sudo make’.)

The code **class_sz** runs in parallel, so you need a **gcc** compiler that is not **clang**.

The previous commands compile both the executable and the python wrapper.
If you do not want to compile the **classy** python module do ‘$ make class’.

For the python module, you need the prerequisites such as **numpy**, **scipy**
and **Cython** installed on your computer.

(class_sz also works on the new mac M1 chips.)

Run the code with most of the power spectra output:

    $ ./class class_sz_test.ini

Run the code with a simple tSZ computation:

    $ ./class class-sz_simple.ini


The  'ini' files are the parameter files. I will be releasing a detailed explanatory file soon.

If any of these two ini files crash, it simply means that the installation was not successful. In this case, please read carefully this readme file and follow the instructions given below. If you are still not able to run these test files, please get in touch.
If nothing appears to solve your installation issues: it is a good idea to try installing the original class code and check that it runs as well as its python wrapper (e.g., the notebook cl_ST.ipynb). If the class code does not run on your system, you should consult the issue page of the class repository and first make sure you solve your issues with the original class code, before moving to class_sz.


Computing SZ and Halo model quantities via the Python wrapper classy_sz
------------------------------


Once class_sz is installed. You can use classy_sz just as you use classy with the normal class code.
You can compute everything classy computes, as well as all the halo model quantities implemented in class_sz.

First, make sure that you have compiled the python wrapper with:

$ make clean

$ make

(Note that the second command must be 'make', and not 'make class' for the python wrappper to be compiled.)

That's it!

Have a look at the notebook class_sz_plots_and_tutorial.ipynb and try to run it. It should output the primary cmb and tsz power spectra.
The notebook is here:

https://github.com/borisbolliet/class_sz/blob/master/notebooks/class_sz_plots_and_tutorial.ipynb


Some tips to run on computer clusters
------------------------------

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

class_sz now requires FFTW3 library, used for the computations of kSZ^2 x LSS power spectra and bispectra.

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
In my case the modif looks like this:

  extra_link_args=['-lgomp','-lgsl','-lgslcblas','**-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/11/**']


Compiler - GCC version
------------------------------

The current gcc version used in the makefile is gcc-11. But this  can be changed easily to any gcc version that is available to you.
There are two modifications:

1) Line 20 of Makefile: CC = gcc-XX (where XX=11 in my case.)

2) Line 12 of python/setup.py: replace 'gcc-11' with, e.g., 'gcc-XX'.



Support
-------

To get support on the class_sz module, feel free to contact me via slack/email (boris.bolliet@gmail.com), or open an issue on the GitHub page.

Acknowledgment
-------

Thanks to  Juan Macias-Perez, Eiichiro Komatsu, Ryu Makiya, Barabara Comis, Julien Lesgourgues, Jens Chluba, Colin Hill, Florian Ruppin, Thejs Brinckmann, Aditya Rotti, Mathieu Remazeilles, David Alonso, Nick Koukoufilippas, Fiona McCarthy, Eunseong Lee, Ola Kusiak, Simone Ferraro, Mat Madhavacheril, Manu Schaan, Shivam Pandey for help, suggestions and/or running tests with **class_sz**.
