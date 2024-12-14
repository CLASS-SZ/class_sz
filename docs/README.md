
# Installation

To install:

```
pip install classy_sz
```

should work on most platforms as we uploaded pre-compiled binaries to [PyPi](https://pypi.org/project/classy-sz/#files) for Linux and MacOS, and Python version >= 3.9.

## setting path for CLASS-SZ data (optional)

By default, the neural nets emulators (~1GB of files) will be installed in your home directory. If you're working on a computing cluster or prefer to store the data elsewhere, you can specify a custom directory.

To specify where you want to store the neural nets data, run the following command in your terminal **before** importing the package:

```bash
export PATH_TO_CLASS_SZ_DATA=/path/to/store/data
mkdir -p $PATH_TO_CLASS_SZ_DATA/class_sz_data_directory
```

This command sets the `PATH_TO_CLASS_SZ_DATA` variable for the current session.

To ensure this variable is set every time you open a terminal, you can add this line to your `~/.bashrc` or `~/.bash_profile` file automatically using the `echo` command.

For `~/.bashrc` (common for most Linux systems), type in your terminal:
```bash
echo -e "\n# Set path for CLASS-SZ data\nexport PATH_TO_CLASS_SZ_DATA=/path/to/store/data" >> ~/.bashrc
echo -e "\n# Create directory for CLASS-SZ data\nmkdir -p \$PATH_TO_CLASS_SZ_DATA/class_sz_data_directory" >> ~/.bashrc
```

To apply the changes immediately:
```bash
source ~/.bashrc
```

(Replace `bashrc` by `bash_profile` if you use macOS.)

Now, every time you open a terminal, the `PATH_TO_CLASS_SZ_DATA` environment variable will automatically be set to your specified directory, ensuring the neural nets emulators are always stored in the correct location.

