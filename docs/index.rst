.. cmbagent documentation master file, created by
   sphinx-quickstart on Sun Sep 15 15:20:33 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

class_sz
====================================

.. image:: https://img.shields.io/github/stars/CLASS-SZ/class_sz?style=social
   :target: https://github.com/CLASS-SZ/class_sz
   :alt: GitHub stars

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: LICENSE
   :alt: License

.. image:: https://readthedocs.org/projects/class_sz/badge/?version=latest
   :target: https://class_sz.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/classy_sz.svg
   :target: https://pypi.org/project/classy_sz/
   :alt: PyPI version


Cosmic Linear Anisotropy Solving System with Machine Learning Accelerated CMB, LSS and Halo Model Observables Computations

.. note:: 

   .. raw:: html

      <div style="border: 1px solid #7cb342; background-color: #e8f5e9; padding: 10px; border-radius: 5px;">
      
   **Class_sz** is compatible with **Jax** and now allows for **automatic differentiation** on some of its output. See
   `here <https://class-sz.readthedocs.io/en/latest/notebooks/classy_szfast_matter_pk_linear.html#Gradients-at-all-k's>`_ 
   for an example on the matter power spectrum. The code can now be used in Hamiltonian Monte Carlo and Simulation-Based 
   Inference pipelines.

   .. raw:: html

      </div>


.. note::
   This software is being actively developed. We welcome contributions from the community.


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   README.md


.. toctree::
   :maxdepth: 2
   :caption: Advanced Installation

   notebooks/classy_szfast_install


.. toctree::
   :maxdepth: 2
   :caption: Input parameters

   notebooks/class_szfast_parameters_and_derived_parameters


   
.. toctree::
   :maxdepth: 2
   :caption: Hubble function

   notebooks/classy_szfast_hz


.. toctree::
   :maxdepth: 2
   :caption: Comoving distance

   notebooks/classy_szfast_daz

.. toctree::
   :maxdepth: 2
   :caption: CMB power spectra 

   notebooks/classy_szfast_cmb_cls
   

.. toctree::
   :maxdepth: 2
   :caption: Linear matter power spectrum

   notebooks/classy_szfast_matter_pk_linear


.. toctree::
   :maxdepth: 2
   :caption: Biased tracers

   notebooks/class_sz_quaia_example


.. toctree::
   :maxdepth: 2
   :caption: Halo mass function

   notebooks/classy_szfast_hmf_and_sigma



.. toctree::
   :maxdepth: 2
   :caption: CMB Lensing x Galaxy lensing

   notebooks/class_sz_cmblensing-x-galaxylensing
   
   
.. toctree::
   :maxdepth: 2
   :caption: Galaxy x Galaxy lensing

   notebooks/class_sz_galaxy-x-galaxylensing
   

.. toctree::
   :maxdepth: 2
   :caption: TSZ power spectrum

   notebooks/class_sz_tszpowerspectrum
   


.. toctree::
   :maxdepth: 2
   :caption: TSZ X Galaxies

   notebooks/class_sz_ngal_tSZ


.. toctree::
   :maxdepth: 2
   :caption: CIB auto

   notebooks/class_sz_cibauto


.. toctree::
   :maxdepth: 2
   :caption: CIB x CMB lensing

   notebooks/class_sz_cibxcmblensing
   

.. toctree::
   :maxdepth: 2
   :caption: kSZ2 X Galaxies

   notebooks/class_sz_kSZ2g


.. toctree::
   :maxdepth: 2
   :caption: SZ cluster cosmology

   notebooks/class_szfast_cluster_counts



   
