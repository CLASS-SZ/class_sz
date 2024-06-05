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

# Function to find the cosmopower-organization directory
find_cosmopower_dir() {
    find "$1" -maxdepth "$DEPTH" -type d -name "cosmopower-organization" -print -quit
}



COSMOPOWER_DIR=$(find_cosmopower_dir "..")

if [ -z "$COSMOPOWER_DIR" ]; then
    COSMOPOWER_DIR=$(find_cosmopower_dir ".")
fi

if [ -n "$COSMOPOWER_DIR" ]; then
    echo "Found cosmopower-organization directory at: $(realpath "$COSMOPOWER_DIR")"
    PATH_TO_COSMOPOWER_ORGANIZATION=$(realpath "$COSMOPOWER_DIR")
else
    echo "--> cosmopower-organization directory not found within $DEPTH levels up or down."
    echo "--> We will install it one level up!"
    cd ..
    mkdir cosmopower-organization
    cd cosmopower-organization
    git clone https://github.com/cosmopower-organization/lcdm.git
    git clone https://github.com/cosmopower-organization/mnu.git
    git clone https://github.com/cosmopower-organization/mnu-3states.git
    git clone https://github.com/cosmopower-organization/ede.git
    git clone https://github.com/cosmopower-organization/neff.git
    git clone https://github.com/cosmopower-organization/wcdm.git
    cd ..
    PATH_TO_COSMOPOWER_ORGANIZATION=$(realpath "cosmopower-organization")
fi


# Write environment setup to file
echo "export PATH_TO_COSMOPOWER_ORGANIZATION=${PATH_TO_COSMOPOWER_ORGANIZATION}" > "${MDIR}/../set_class_sz_env.sh"
echo "PATH_TO_COSMOPOWER_ORGANIZATION is set to ${PATH_TO_COSMOPOWER_ORGANIZATION}"

# Export the variable
export PATH_TO_COSMOPOWER_ORGANIZATION="${PATH_TO_COSMOPOWER_ORGANIZATION}"
