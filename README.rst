==============================================
CLASS_SZ
==============================================
 Cosmic Linear Anisotropy Solving System

 with thermal Sunyaev Zeldovich power spectrum computation

 and other halo model observables

**The SZ module is initially based on Eiichiro Komatsu’s fortran code SZFAST.**

(See http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/)


The SZ module is included in the file **source/szpowerspectrum.c**
and its dependencies.

In addition to SZ power spectrum, class_sz can compute cross and auto power spectra for other tracers
in the halo model (currently being developped: kSZ, galaxy, ISW, lensing and CIB).

The code is currently in development, don't hesitate to reach out if anything is unclear due to lack of comments and indications.

**The code CLASS_SZ is an extension of the CLASS code.**

For download and information on **class**, see http://class-code.net and https://github.com/lesgourg/class_public


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

Run the code with SZ power spectrum computation

    $ ./class class-sz_simple.ini


The explanatory files are reference input files, containing and
explaining the use of all possible input parameters.


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

Note that these prescription are system dependent: you may not need them if your path and environment variables are such that gsl and its libraries are well linked.

MacOS problem with OpenMP
------------------------------

To run the code in parallel, you may run into a problem on a mac. The solution is provided here:

https://github.com/lesgourg/class_public/issues/208

Essentially, you need to edit a line in python/setup.py such as the code knows about the mpi libraries to be used with your compiler (gcc-10 in the example below).
In my case the modif looks like this:

  extra_link_args=['-lgomp','-lgsl','-lgslcblas','**-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/10/**']


Compiler - GCC version
------------------------------

The current gcc version used in the makefile is gcc-10. But this  can be changed easily to any gcc version that is available to you.
There is two modifications:

1) Line 20 of Makefile: CC = gcc-XX (where XX=10 in my case.)

2) Line 12 of python/setup.py: replace 'gcc-10' with, e.g., 'gcc-XX'.


Computing SZ and Halo model quantities via the Python wrapper classy
------------------------------


Once class_sz is installed. You can use classy just as you use it with the normal class code.

First, make sure that you have compiled the python wrapper with:

$ make clean

$ make

(Note that the second command must be 'make', and not 'make class' for the python wrappper to be compiled.)

Then, the only extra-requirement compared to the standard class is that the following path needs to be set in the dictionary of parameters passed to classy.

  path_to_class : '/path/to/class_sz'

That's it!

Have a look at the notebook cl_ST_sz.ipynb and try to run it. It should output the primary cmb and tsz power spectra.


Computing and Plotting results
------------------------------

There is a script that enables to plot tSZ spectra easily, e.g.:


    /path/to/python /path/to/class_sz/sz_auxiliary_files/run_scripts/tSZ_varying_params.py -param_name M1SZ -min 1e13 -max 1e15 -N 3 -spacing log -output 'tSZ_1h,tSZ_Trispectrum'  -show_legend yes -show_error_bars yes -compute_scaling_with_param no -save_tsz_ps no -plot_ref_data no -print_rel_diff no


One just needs to set the path to class_sz properly at the beginning of tSZ_varying_params.py, and it should run.



Support
-------

To get support on the SZ module, feel free to contact me via slack/email (boris.bolliet@gmail.com), or open an issue on the GitHub page.

Acknowledgment
-------

Thanks to  Juan Macias-Perez, Eiichiro Komatsu, Ryu Makiya, Barabara Comis, Julien Lesgourgues, Jens Chluba, Colin Hill, Florian Ruppin, Thejs Brinckmann, Aditya Rotti, Mathieu Remazeilles, David Alonso, Nick Koukoufilippas, Fiona McCarthy for help, suggestions and/or running tests with **class_sz**.
