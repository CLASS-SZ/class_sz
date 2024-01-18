#  halomodel_cib_tsz_cibxtsz


This code has been been adapted to run specifically with class_sz.

For the original code see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz

To run the benchmarking against class_sz do:

$ cd sz_auxiliary_files/maniyar_cib_benchmark/

$ python driver_cell.py


Computes the CIB, tSZ, and CIB-tSZ correlation power spectrum using a newly developed halo model

* This code is based on a newly developed halo model for the CIB and CIB-tSZ correlation with just four physical model parameters.
* These models can be used in the CMB data analysis to account for the CIB, tSZ, and CIBxtSZ foregrounds instead of just fitting them with a power law or with some templates.
* hmf_unfw_bias.py code is used to calculate the halo mass function, Fourier transform of the NFW profile, and the halo bias. The code to calculate the halo mass function has been adopted from the 'hmfcalc' code available online: https://github.com/steven-murray/hmf.

Clone the repository, then compute the power spectra with:
```
python driver_cell.py
```
or look at ```driver_cell.ipynb```.

Some things to be careful about:

* Matter power spectra for redshifts from 'redshifts.txt' are stored in files and are directly read.
However, their numbering is backwards i.e. _1 file has highest redshift and _210 has lowest. An easier
option might be to create your own matter power spectra from camb and read them in the code directly.
* For Planck and Herschel, the SEDs filtered with their respective bandpasses are provided. For other
experiments, the code reads raw SEDs and calculates the power spectra for given frequencies. You can
of course replace the raw SEDs with filtered SEDs of your instrument of choice if they are available.
* Just like the matter power spectra, if you want to use your own code/ some other code to generate
the halo mass function/Fourier tranform of the NFW profile/halo bias, that is possible as well.
* All the changes proposed above will have to be done in the ```input_var.py``` file.

Hope you find this code useful! If you use this code in a publication, please cite https://arxiv.org/abs/2006.16329.
Do not hesitate to contact me with any questions: abhimaniyar@gmail.com
