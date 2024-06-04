import os
import sys
import subprocess as sbp
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import numpy as np
import site

# Get the site-packages path
path_install = site.getusersitepackages() if os.path.exists(site.getusersitepackages()) else site.getsitepackages()[0]
path_install = os.path.join(path_install, "class-sz")

# Helper function to include additional package files
def package_files(directory, exclude_dirs):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        # Exclude directories
        directories[:] = [d for d in directories if os.path.join(path, d) not in exclude_dirs]
        for filename in filenames:
            paths.append(os.path.relpath(os.path.join(path, filename), directory))
    return paths

# Define directories to exclude
exclude_dirs = [
    os.path.join('class-sz', 'output'),
    os.path.join('class-sz', 'ORTHOGONAL_NG_matrices'),
    os.path.join('class-sz', 'pt_matrices')
]

# Gather all package files, excluding specific directories
pck_files = package_files('class-sz', exclude_dirs)

# Get the GCC compiler path
GCCPATH_STRING = sbp.Popen(['gcc', '-print-libgcc-file-name'], stdout=sbp.PIPE).communicate()[0]
GCCPATH = os.path.normpath(os.path.dirname(GCCPATH_STRING)).decode()

# Determine the libraries to link against
liblist = ["class"]
MVEC_STRING = sbp.Popen(['gcc', '-lmvec'], stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec", "m"]

# Define absolute paths
root_folder = "class-sz"
include_folder = os.path.join(root_folder, "include")
classy_sz_folder = os.path.join(root_folder, "python")

# Define the Cython extension
classy_sz_ext = Extension(
    "classy_sz",
    [os.path.join(classy_sz_folder, "classy.pyx")],
    include_dirs=[np.get_include(), include_folder],
    libraries=liblist,
    library_dirs=[root_folder, GCCPATH],
    extra_link_args=['-lgomp', '-lgsl', '-lfftw3', '-lgslcblas']
)
classy_sz_ext.cython_directives = {'language_level': "3" if sys.version_info.major >= 3 else "2"}

# Custom build_ext class
class ClassyBuildExt(build_ext):
    def run(self):
        run_env = dict(CLASSDIR=path_install, **os.environ.copy())
        
        # Run the script to select the correct Makefile and build
        sbp.run(["chmod", "+x", "select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)
        sbp.run(["./select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)
        
        # Proceed with the standard build_ext behavior
        build_ext.run(self)

# Long description for the package
long_description = "The official repository for the classy_sz code ('https://github.com/CLASS-SZ')."

# Setup function
setup(
    name='classy_sz',
    version="0.0.7",
    author="Boris Bolliet, Ola Kusiak",
    author_email="bb667@cam.ac.uk, akk2175@columbia.edu",
    description='CLASS-SZ in Python',
    long_description=long_description,
    url='https://github.com/CLASS-SZ',
    cmdclass={'build_ext': ClassyBuildExt},
    ext_modules=[classy_sz_ext],
    packages=find_packages(where='class-sz/python', exclude=['class-sz.output', 
                                                            'class-sz.ORTHOGONAL_NG_matrices',
                                                            'class-sz.pt_matrices']),
    package_dir={'': 'class-sz/python'},
    package_data={'': pck_files},
    include_package_data=True,
    install_requires=["classy_szfast>=0.0.7"],
)