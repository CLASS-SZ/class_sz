import os
import sys
import subprocess as sbp
from setuptools import setup, find_packages, Extension
# from setuptools.command.install import install
from Cython.Distutils import build_ext
import numpy as np
import site
from setuptools.command.install import install

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
exclude_dirs = []

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
        sbp.Popen("make libclass.a -j",shell=True,cwd=os.path.join(os.getcwd(),"class-sz"),env=run_env).wait()
        sbp.run(["chmod", "+x", "download_emulators.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)
        sbp.Popen("./download_emulators.sh",shell=True, cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env).wait()
        
        # Proceed with the standard build_ext behavior
        build_ext.run(self)

class CustomInstallCommand(install):
    """Customized setuptools install command to run additional steps after installation."""
    def run(self):
        # Run the standard install
        install.run(self)

        # Check if the environment variable for data directory is set
        data_dir = os.environ.get('CLASS_SZ_DATA_DIR')

        if data_dir:
            # Check if the data directory exists and is not empty
            if os.path.exists(data_dir) and os.listdir(data_dir):
                print(f"Data directory {data_dir} already exists and is not empty. Skipping data download.")
            else:
                # Create the data directory if it doesn't exist
                os.makedirs(data_dir, exist_ok=True)
                # Define the GitHub repository
                for cosmo_model in ['lcdm', 'mnu', 'mnu-3states', 'ede', 'neff', 'wcdm']:
                    print(f"Downloading data for {cosmo_model}")
                    github_repo = f"https://github.com/cosmopower-organization/{cosmo_model}.git"
                    
                    # Clone the GitHub repository to a temporary directory
                    temp_dir = os.path.join(os.getcwd(), "temp_repo")
                    if not os.path.exists(temp_dir):
                        os.makedirs(temp_dir)
                    subprocess.check_call(['git', 'clone', github_repo, temp_dir])
                    
                    # Move the required data to the specified directory
                    # Assuming the data is located in 'data' directory in the cloned repo
                    data_source_dir = os.path.join(temp_dir, 'data')
                    if os.path.exists(data_source_dir):
                        for item in os.listdir(data_source_dir):
                            source = os.path.join(data_source_dir, item)
                            destination = os.path.join(data_dir, item)
                            if os.path.isdir(source):
                                subprocess.check_call(['cp', '-r', source, destination])
                            else:
                                subprocess.check_call(['cp', source, destination])
                    
                    # Clean up the temporary directory
                    subprocess.check_call(['rm', '-rf', temp_dir])
                
                print(f"Data downloaded to {data_dir}")
        else:
            print("CLASS_SZ_DATA_DIR environment variable is not set. Skipping data download.")

long_description = "See ('https://github.com/CLASS-SZ')."

# Setup function
setup(
    name='classy_sz',
    version="0.1.1",
    author="Boris Bolliet, Ola Kusiak",
    author_email="bb667@cam.ac.uk, akk2175@columbia.edu",
    description='CLASS-SZ in Python',
    long_description=long_description,
    url='https://github.com/CLASS-SZ',
    cmdclass={'build_ext': ClassyBuildExt, 'install': CustomInstallCommand},
    ext_modules=[classy_sz_ext],
    packages=find_packages(where='class-sz/python'),
    package_dir={'': 'class-sz/python'},
    package_data={'': pck_files},
    include_package_data=True,
)