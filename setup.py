import os
import sys
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import numpy as np
import site
import subprocess as sbp

# Get the site-packages path
path_install = site.getusersitepackages() if os.path.exists(site.getusersitepackages()) else site.getsitepackages()[0]
path_install = os.path.join(path_install, "class-sz")

# Helper function to include additional package files
def package_files(directory, exclude_dirs):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        directories[:] = [d for d in directories if os.path.join(path, d) not in exclude_dirs]
        for filename in filenames:
            paths.append(os.path.relpath(os.path.join(path, filename), directory))
    return paths

exclude_dirs = ['class-sz/output', 'class-sz/test', 'class-sz/class_sz_auxiliary_files/excludes']
pck_files = package_files('class-sz', exclude_dirs)

# Get the GCC compiler path
GCCPATH_STRING = sbp.Popen(['gcc', '-print-libgcc-file-name'], stdout=sbp.PIPE).communicate()[0]
GCCPATH = os.path.normpath(os.path.dirname(GCCPATH_STRING)).decode()

# Determine libraries to link
liblist = ["class"]
MVEC_STRING = sbp.Popen(['gcc', '-lmvec'], stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec", "m"]

# Check the platform to determine the correct OpenMP library
openmp_flag = '-lgomp' if sys.platform.startswith('linux') else '-lomp'


# Define Cython extension
classy_sz_ext = Extension(
    "classy_sz",
    [os.path.join("class-sz/python", "classy.pyx")],
    include_dirs=[np.get_include(), os.path.join("class-sz", "include")],
    libraries=liblist,
    library_dirs=["class-sz", GCCPATH],
    extra_link_args=[openmp_flag, '-lgsl', '-lfftw3', '-lgslcblas']
)

classy_sz_ext.cython_directives = {'language_level': "3" if sys.version_info.major >= 3 else "2"}

# Custom build_ext class to handle building `libclass.a`
class ClassyBuildExt(build_ext):
    def run(self):
        run_env = dict(CLASSDIR=path_install, **os.environ.copy())
        # Print the current working directory for debugging purposes
        print("Current working directory:", os.getcwd())

        # Change permissions on the select_makefile.sh script
        print("Running chmod to make select_makefile.sh executable")
        sbp.run(["chmod", "+x", "select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)

        # Run the script to select the correct Makefile
        print("Running the script to select the correct Makefile")
        sbp.run(["./select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)

        print("Building the library") ## prints wont print unless -vvv is used
        result = sbp.run("make libclass.a -j", shell=True, cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env)
        if result.returncode != 0:
            raise RuntimeError("Building libclass.a failed")

        build_ext.run(self)

print(pck_files)
# Setup function
setup(
    name='classy_sz',
    version="0.1.45",
    author="Boris Bolliet, Ola Kusiak",
    author_email="bb667@cam.ac.uk, akk2175@columbia.edu",
    description='CLASS-SZ in Python',
    long_description="See ('https://github.com/CLASS-SZ').",
    url='https://github.com/CLASS-SZ',
    cmdclass={'build_ext': ClassyBuildExt},
    ext_modules=[classy_sz_ext],
    packages=find_packages(where='class-sz/python'),
    package_dir={'': 'class-sz/python'},
    package_data={'class-sz': pck_files},
    include_package_data=True,
)
