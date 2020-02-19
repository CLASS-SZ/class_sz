First steps towards an SZ power spectrum likelihood for SO

To run the SZ power spectrum likelihood, Cobaya needs a small edit.

Open the file 'cobaya/theories/classy/classy.py' of the Cobaya directory
Paste the function 'get_Cl_sz(self)' from class_sz/sz_auxiliary_files/cobaya_class_sz_likelihoods/cobaya_classy_file/classy_modif.dat
into 'cobaya/theories/classy/classy.py' [see classy_modif.dat]

Then, with the terminal move to the cobaya directory and run:

$ python setup.py install

Then, to make sure cobaya will use class_sz rather than class, just 'cd' to class_sz directory and run:

$ make clean
$ make

So that now the 'classy' used by cobaya is actually class_sz !

The SZ power spectrum likelihood code is called 'so_sz_y.py' and is located in 'class_sz/sz_auxiliary_files/cobaya_class_sz_likelihoods/likelihoods'

To run the likelihood code you can use 'run_so_sz.py' located in the same directory, by simply doing:

$ python run_so_sz.py
