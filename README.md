
# CLASS_SZ

Cosmic Linear Anisotropy Solving System with Machine Learning Accelerated and Accurate CMB, LSS, and Halo Model Observables Computations

To install/run CLASS_SZ on the fly, check the Colab notebook:

[Colab Notebook](https://colab.research.google.com/drive/1AULgG4ZLLG1YXRI86L54-hpjWyl1X-8c?usp=sharing).

Further information on downloading and installing the code is given hereafter for Mac, M1/2 Mac, and Linux.

The tutorial notebooks can be found at:

[CLASS_SZ Notebooks](https://github.com/CLASS-SZ/notebooks)

These notebooks along with the paper ([Bolliet et al 2023](https://arxiv.org/abs/2310.18482)) constitute the current documentation.

CLASS_SZ is as fast as it gets, with full parallelization, implementation of high-accuracy neural network emulators, and Fast Fourier Transforms.

Since it is based on Lesgourgues's class code, the halo model and LSS calculations (essentially based on distances and matter clustering) are always consistent with the cosmological model computed by CLASS.

CLASS_SZ is an extension of Julien Lesgourgues's [CLASS](https://github.com/lesgourg/class_public) code.

CLASS_SZ is initially based on Eiichiro Komatsu’s Fortran code [SZFAST](http://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumks/).

CLASS_SZ functionalities are located in the files:

- **source/class_sz.c** for the main CLASS_SZ functions,
- **tools/class_sz_tools.c** for useful routines,
- **source/class_sz_clustercounts.c** for tSZ cluster counts. Since March 2024, CLASS_SZ cluster counts calculations are superseded by [cosmocnc](https://github.com/inigozubeldia/cosmocnc) ([Zubeldia & Bolliet 2024](https://arxiv.org/abs/2403.09589)).

CLASS_SZ's outputs are regularly cross-checked with other CMBxLSS codes, such as:

- [cosmocnc](https://github.com/inigozubeldia/cosmocnc)
- [hmvec](https://github.com/simonsobs/hmvec/tree/master/hmvec)

## Installation

To compile the code, you need a C compiler, a Fortran compiler, and FFTW. To get these, do the following:

### MacOS

Install Xcode and its command-line tools:

```sh
xcode-select --install
```

Install Homebrew if you don't have it:

```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Install dependencies:

```sh
brew install gcc fftw
```

### Linux

Install dependencies:

```sh
sudo apt-get update
sudo apt-get install build-essential gfortran libfftw3-dev
```

## Compilation

To compile CLASS_SZ, simply run:

```sh
make
```

## Running CLASS_SZ

To run CLASS_SZ, you need an input parameter file. An example is provided in `explanatory.ini`.

Run the code as follows:

```sh
./class explanatory.ini
```

## Citing CLASS_SZ

If you use CLASS_SZ in your research, please cite the following paper:

[Bolliet et al 2023](https://arxiv.org/abs/2310.18482)

## Contributing

Contributions are welcome. Please submit pull requests or report issues on GitHub.

## License

CLASS_SZ is distributed under the terms of the MIT license.
