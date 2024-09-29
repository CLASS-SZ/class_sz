
# Change to the class-sz directory
cd class-sz/

# Clean previous builds
make clean

# Go back to the root directory
cd ..

# Remove old build files
rm -rf build dist *.egg-info
rm -rf class-sz/python/classy_sz.*so
rm -rf class-sz/python/classy_sz.*egg-info

# Check if Docker is running (this is relevant for Linux builds)
docker --version

# Set environment variables and build wheels for Windows
CIBW_BUILD="cp310-*" \
CIBW_SKIP="cp27-* cp34-* cp35-* cp36-* cp37-* cp38-* pp* *-win32" \
CIBW_ARCHS="x86" \
CIBW_BEFORE_BUILD="curl -O https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz && \
  tar -xzf gsl-2.7.1.tar.gz && \
  cd gsl-2.7.1 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && \
  curl -O http://www.fftw.org/fftw-3.3.10.tar.gz && \
  tar -xzf fftw-3.3.10.tar.gz && \
  cd fftw-3.3.10 && \
  ./configure --prefix=/usr/local --enable-shared && \
  make && \
  make install" \

CIBW_PLATFORM=windows

python -m cibuildwheel --output-dir dist

