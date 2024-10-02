#!/bin/bash

# Get current directory
MDIR=$(pwd)
WRKDIR="${MDIR}/build"
OUTDIR="${MDIR}/output"

# Ensure directories exist
if [ ! -e "${WRKDIR}" ]; then
    mkdir -p "${WRKDIR}/lib"
fi

if [ ! -e "${OUTDIR}" ]; then
    mkdir "${OUTDIR}"
fi

touch build/.base

# Depth for searching directories
DEPTH=3

# Function to find the class_sz_data_directory directory
find_class_sz_data_dir() {
    find "$1" -maxdepth "$DEPTH" -type d -name "class_sz_data_directory" -print -quit
}



CLASS_SZ_DATA_DIR=$(find_class_sz_data_dir "..")

if [ -z "$CLASS_SZ_DATA_DIR" ]; then
    CLASS_SZ_DATA_DIR=$(find_class_sz_data_dir ".")
fi

if [ -n "$CLASS_SZ_DATA_DIR" ]; then
    echo "Found class_sz_data_directory directory at: $(realpath "$CLASS_SZ_DATA_DIR")"
    PATH_TO_CLASS_SZ_DATA=$(realpath "$CLASS_SZ_DATA_DIR")
else
    echo "--> class_sz_data_directory directory not found within $DEPTH levels up or down."
    echo "--> We will install it one level up!"
    cd ..
    mkdir class_sz_data_directory
    cd class_sz_data_directory
    git clone https://github.com/cosmopower-organization/lcdm.git
    git clone https://github.com/cosmopower-organization/mnu.git
    git clone https://github.com/cosmopower-organization/mnu-3states.git
    git clone https://github.com/cosmopower-organization/ede.git
    git clone https://github.com/cosmopower-organization/neff.git
    git clone https://github.com/cosmopower-organization/wcdm.git
    cd ..
    PATH_TO_CLASS_SZ_DATA=$(realpath "class_sz_data_directory")
fi


# Write environment setup to file
echo "export PATH_TO_CLASS_SZ_DATA=${PATH_TO_CLASS_SZ_DATA}" > "${MDIR}/../set_class_sz_env.sh"
echo "PATH_TO_CLASS_SZ_DATA is set to ${PATH_TO_CLASS_SZ_DATA}"

# Export the variable
export PATH_TO_CLASS_SZ_DATA="${PATH_TO_CLASS_SZ_DATA}"