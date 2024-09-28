#!/bin/bash

# Ensure to manually update the version number in setup.py before running this script

# Function to increment version
increment_version() {
  current_version=$(grep -Eo 'version="[^"]*' setup.py | cut -d'"' -f2)
  IFS='.' read -r -a version_parts <<< "$current_version"

  major=${version_parts[0]}
  minor=${version_parts[1]}
  patch=${version_parts[2]}

  patch=$((patch + 1))
  if [ $patch -ge 100 ]; then
    patch=0
    minor=$((minor + 1))
  fi
  if [ $minor -ge 10 ]; then
    minor=0
    major=$((major + 1))
  fi

  new_version="$major.$minor.$patch"

  # Update the setup.py with the new version
  sed -i '' "s/version=\"$current_version\"/version=\"$new_version\"/" setup.py
}

# Increment the version
# increment_version

cd class-sz/

make clean

cd ..

# Remove old build files
rm -rf build dist *.egg-info

rm -rf class-sz/python/classy_sz.*so
rm -rf class-sz/python/classy_sz.*egg-info

# check docker is running
docker --version

# build wheels for linux
# CIBW_BEFORE_BUILD="apt-get install -y libgsl-dev libfftw3-dev" CIBW_PLATFORM=linux cibuildwheel --output-dir dist
# CIBW_BUILD="cp39-manylinux_aarch64" CIBW_ARCHS="aarch64" CIBW_BEFORE_BUILD="yum install -y gsl-devel fftw-devel" CIBW_PLATFORM=linux cibuildwheel --output-dir dist


# 3. Windows Wheels
CIBW_BUILD="cp39-win_amd64 cp39-win32" \
CIBW_PLATFORM=windows cibuildwheel --output-dir dist

# # build the wheel for mac
# # Build the package
# python setup.py sdist bdist_wheel

# # Upload the package
# twine upload dist/*

# # Clean up build directories
# rm -rf build dist *.egg-info