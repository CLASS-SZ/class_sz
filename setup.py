import os
import sys
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import numpy as np
import site
import subprocess as sbp
import platform
import shutil


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

# Check for GCC on Windows (via MinGW or WSL)
def get_gcc_path():
    try:
        gcc_output = sbp.check_output(['gcc', '-print-libgcc-file-name']).strip()
        gcc_path = os.path.normpath(os.path.dirname(gcc_output))
        return gcc_path.decode() if isinstance(gcc_path, bytes) else gcc_path
    except Exception as e:
        print(f"Error finding GCC: {e}")
        return None

GCCPATH = get_gcc_path()

# Libraries to link
liblist = ["class"]
try:
    MVEC_STRING = sbp.check_output(['gcc', '-lmvec'], stderr=sbp.PIPE)
    if b"mvec" not in MVEC_STRING:
        liblist += ["mvec", "m"]
except Exception:
    # GCC or mvec might not exist on Windows; handle it gracefully
    liblist += ["m"]

# Check the platform to determine the correct OpenMP library
openmp_flag = '-lgomp' if sys.platform.startswith('linux') else '-lomp'

# Define Cython extension
classy_sz_ext = Extension(
    "classy_sz",
    [os.path.join("class-sz/python", "classy.pyx")],
    include_dirs=[np.get_include(), os.path.join("class-sz", "include")],
    libraries=liblist,
    library_dirs=["class-sz"] + ([GCCPATH] if GCCPATH else []),
    extra_link_args=[openmp_flag, '-lgsl', '-lfftw3', '-lgslcblas']
)

classy_sz_ext.cython_directives = {'language_level': "3" if sys.version_info.major >= 3 else "2"}

# # Custom build_ext class to handle building `libclass.a`
# class ClassyBuildExt(build_ext):
#     def run(self):
#         run_env = dict(CLASSDIR=path_install, **os.environ.copy())
#         print("Current working directory:", os.getcwd())

#         # Skip chmod on Windows since it's unnecessary
#         if not sys.platform.startswith('win'):
#             print("Running chmod to make select_makefile.sh executable")
#             sbp.run(["chmod", "+x", "select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)

#         # Run the script to select the correct Makefile
#         print("Running the script to select the correct Makefile")
#         if os.path.exists(os.path.join(os.getcwd(), "class-sz", "select_makefile.sh")):
#             print("Running the script to select the correct Makefile")
#             sbp.run(["./select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)
#         else:
#             raise FileNotFoundError("select_makefile.sh not found in class-sz directory")

#         sbp.run(["./select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)
#         # sbp.run(["bash", "./select_makefile.sh"], cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env, check=True)

        
#         # Build the library
#         print("Building the library")
#         result = sbp.run("make libclass.a -j", shell=True, cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env)
#         if result.returncode != 0:
#             raise RuntimeError("Building libclass.a failed")

#         build_ext.run(self)



# Custom build_ext class to handle building `libclass.a`
class ClassyBuildExt(build_ext):
    def run(self):
        run_env = dict(CLASSDIR=path_install, **os.environ.copy())
        print("Current working directory:", os.getcwd())

        def select_makefile():
            system = platform.system()
            machine = platform.machine()

            makefile_dir = os.path.join(os.getcwd(), "class-sz")  # Ensure correct directory
            print(f"System: {system}")
            print(f"Machine: {machine}")
            print(f"Makefile directory: {makefile_dir}")

            # Select the appropriate Makefile based on OS and architecture
            if system == "Darwin" and machine == "arm64":
                print("Using M1 Makefile")
                shutil.copy(os.path.join(makefile_dir, "Makefile_m1"), os.path.join(makefile_dir, "Makefile"))
            elif system == "Linux":
                print("Using Linux Makefile")
                shutil.copy(os.path.join(makefile_dir, "Makefile_linux"), os.path.join(makefile_dir, "Makefile"))
                print("Environment Variables:")
                print(f"C_INCLUDE_PATH: {os.environ.get('C_INCLUDE_PATH', '')}")
                print(f"LD_LIBRARY_PATH: {os.environ.get('LD_LIBRARY_PATH', '')}")
                print(f"LIBRARY_PATH: {os.environ.get('LIBRARY_PATH', '')}")
                print(f"CFLAGS: {os.environ.get('CFLAGS', '')}")
                print(f"LDFLAGS: {os.environ.get('LDFLAGS', '')}")
            elif system == "Windows":
                print("Using Windows Makefile")
                shutil.copy(os.path.join(makefile_dir, "Makefile_windows"), os.path.join(makefile_dir, "Makefile"))
                print("Environment Variables:")
                print(f"C_INCLUDE_PATH: {os.environ.get('C_INCLUDE_PATH', '')}")
                print(f"LIBRARY_PATH: {os.environ.get('LIBRARY_PATH', '')}")
                print(f"CFLAGS: {os.environ.get('CFLAGS', '')}")
                print(f"LDFLAGS: {os.environ.get('LDFLAGS', '')}")
            else:
                print(f"OS not supported yet: {system}. See README.md for pre-M1 mac support or get in touch.")
                shutil.copy(os.path.join(makefile_dir, "Makefile_m1"), os.path.join(makefile_dir, "Makefile"))

        # Call the Python function to select the correct Makefile
        select_makefile()

        # Build the library using `make` command
        print("Building the library")
        if platform.system() == "Windows":
            # On Windows, ensure MinGW or MSYS2 is installed for `make`
            # Modify the command if necessary to point to the correct `make`
            make_command = "mingw32-make" if shutil.which("mingw32-make") else "make"
        else:
            make_command = "make"

        try:
            result = sbp.run(f"{make_command} libclass.a -j", shell=True, cwd=os.path.join(os.getcwd(), "class-sz"), env=run_env)
            if result.returncode != 0:
                raise RuntimeError("Building libclass.a failed")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            raise RuntimeError(f"Make command not found. Ensure you have `make` installed on your system: {e}")

        # Call the original build_ext run method
        build_ext.run(self)


print(pck_files)
# Setup function
setup(
    name='classy_sz',
    version="0.1.70.post17",
    author="Boris Bolliet, Ola Kusiak, Fiona McCarthy, Alina Sabyr, Kristen Surrao",
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
