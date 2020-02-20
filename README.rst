==============================================
CLASS_SZ
==============================================
 Cosmic Linear Anisotropy Solving System

 With Thermal Sunyaev Zeldovich Power Spectrum Computation


**The SZ module is based on Eiichiro Komatsu’s fortran code SZFAST.**

(See http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/)


The SZ module is included in the file **source/szpowerspectrum.c**
and its dependencies.


**The code CLASS_SZ is an extension of the CLASS code.**

For download and information on **class**, see http://class-code.net and https://github.com/lesgourg/class_public


(README file adapted from the README_CLASS.rst file.)


Downloading the code
--------------

Go to https://github.com/borisbolliet/class_sz_public


Using the code
--------------

You can use **class_sz** freely, provided that in your publications you cite
at least the papers:

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
If you do not want the **classy** python module do ‘$ make class’.

For the python module, you need the prerequisites such as **numpy**, **scipy**
and **Cython** installed on your computer.

Run the code with SZ power spectrum and cluster counts computation

    $ ./class explanatory-sz.ini


The explanatory files are reference input files, containing and
explaning the use of all possible input parameters.

The computation of the tSZ angular power spectrum is stable with masses up to 1e16 Msun/h.


Python Wrapper and Jupyter Notebooks
--------------------------------------

TBD


Support
-------

To get support on the SZ module, feel free to contact me via slack/email
