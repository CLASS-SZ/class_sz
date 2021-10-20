==============================================
CLASS_SZ
==============================================
 Cosmic Linear Anisotropy Solving System

 with fast and accurate halo model computations



In addition to SZ power spectrum, class_sz can compute cross and auto power spectra for other tracers
in the halo model (kSZ, galaxy, ISW, lensing and CIB).

It has several mass functions implemented, with several possible halo mass definitions and concentration-mass
relations. For galaxy clustering and lensing, class_sz has an implementation of HOD based on the one used by
the DES collaboration.

The code is close to be as fast as it can get, with full parallelization.

Since it is based on Lesgourgues's class code, the halo model (essentially based on distances and
matter clustering) is always consistent with the cosmological model.



**Take a look at the notebook to see what class_sz can do:**

https://github.com/borisbolliet/class_sz/blob/master/notebooks/class_sz_plots_and_tutorial.ipynb

The code is currently in development, don't hesitate to reach out if anything is unclear due to lack of comments and indications, or if it crashes unexpectedly.

CLASS_SZ is an extension of Julien Lesgourgues's CLASS code.

For download and information on CLASS, see http://class-code.net and https://github.com/lesgourg/class_public

CLASS_SZ is initially based on Eiichiro Komatsu’s fortran code SZFAST.

(See http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/)

CLASS_SZ modules are located in the files **source/class_sz.c** and **source/szclustercount.c**.


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

References:

`Including massive neutrinos in thermal Sunyaev Zeldovich power spectrum and cluster counts analyses (Bolliet, Brinckmann, Chluba, Lesgourgues, 2020) <https://arxiv.org/abs/1906.10359>`_.

`Dark Energy from the Thermal Sunyaev Zeldovich Power Spectrum (Bolliet, Comis, Komatsu, Macias-Perez, 2017)
<https://arxiv.org/abs/1712.00788>`_.

`CLASS II: Approximation schemes (Blas, Lesgourgues, Tram, 2011)
<http://arxiv.org/abs/1104.2933>`_.

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

Run the code with most of the power spectra output:

    $ ./class class-sz_test.ini

Run the code with a simple tSZ computation:

    $ ./class class-sz_simple.ini


The  'ini' files are the parameter file.
I will be releasing a detailed explanatory file soon.


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



GSL library
------------------------------


New version of class_sz requires gsl (for the integration routines).
One may need to edit the **Makefile** adding the include path for gsl libraries, e.g.,:


    INCLUDES = -I../include -I/usr/local/include/ **-I/path_to_gsl/gsl-2.6/include/**

    class: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS) $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -lm **-L/path_to_gsl/gsl-2.6/lib/ -lgsl -lgslcblas**

For the python wrapper, one also may need to add the absolute path to gsl libraries, e.g.,:

in **class_sz/python/setup.py**:

    classy_ext = Extension("classy", [os.path.join(classy_folder, "classy.pyx")], include_dirs=[nm.get_include(), include_folder, '**/path/to/gsl-2.6/include**'], libraries=liblist,library_dirs=[root_folder, GCCPATH],extra_link_args=['-lgomp','**-L/path_to_gsl/gsl-2.6/lib/**','**-lgsl**','**-lgslcblas**'])



When running, the gsl library also need to be included in the environment variables, i.e., one may
need to do:

    $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path_to_gsl/gsl-2.6/lib

    $ export LD_LIBRARY_PATH

Note that these prescriptions are system dependent: you may not need them if your path and environment variables are such that gsl and its libraries are well linked.

FFTLog library
------------------------------

class_sz now requires FFTW3 library, used for the computations of kSZ^2 x LSS power spectra and bispectra.


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

Thanks to  Juan Macias-Perez, Eiichiro Komatsu, Ryu Makiya, Barabara Comis, Julien Lesgourgues, Jens Chluba, Colin Hill, Florian Ruppin, Thejs Brinckmann, Aditya Rotti, Mathieu Remazeilles, David Alonso, Nick Koukoufilippas, Fiona McCarthy, Eunseong Lee, Ola Kusiak, Simone Ferraro, Mat Madhavacheril, Manu Schaan, for help, suggestions and/or running tests with **class_sz**.
