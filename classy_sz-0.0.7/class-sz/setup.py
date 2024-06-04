from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        name="classy_sz",
        sources=[
            "python/classy.pyx",
            "source/class_sz.c",
            "source/class_sz_clustercounts.c",
            "tools/class_sz_tools.c",
            # Add other source files as necessary
        ],
        include_dirs=[numpy.get_include(), "include"],
        libraries=["gsl", "gslcblas", "fftw3"],
        library_dirs=["/usr/local/lib"],  # Adjust this as necessary
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
]

setup(
    name="class_sz",
    version="0.1",
    description="Cosmic Linear Anisotropy Solving System with Machine Learning Accelerated and Accurate CMB, LSS, and Halo Model Observables Computations",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author="Boris Bolliet",
    url="https://github.com/CLASS-SZ/class_sz",
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    zip_safe=False,
    install_requires=[
        "cython",
        "numpy",
        "scipy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)