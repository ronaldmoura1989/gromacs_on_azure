#!/bin/bash
export DEBIAN_FRONTEND=noninteractive

# Install dependencies and GCC versions
apt-get update -y
apt-get install -y software-properties-common
add-apt-repository universe -y
apt-get update -y
apt-get install -y build-essential cmake libfftw3-dev
apt-get install -y gcc-9 g++-9

# Prepare build directory
WORKDIR="/mnt/gromacs_build"
mkdir -p $WORKDIR
cd $WORKDIR

# Download GROMACS 2023.2
wget -q https://ftp.gromacs.org/gromacs/gromacs-2023.2.tar.gz
tar -zxvf gromacs-2023.2.tar.gz
cd gromacs-2023.2

mkdir build
cd build

# Compile with GCC-9 for CPU (No GPU on D-series)
# Removed -DGMX_GPU=CUDA
cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DREGRESSIONTEST_DOWNLOAD=ON \
  -DCMAKE_C_COMPILER=gcc-9 \
  -DCMAKE_CXX_COMPILER=g++-9 \
  -DGMX_GPU=OFF

make -j64
make install

# Configure environment for the target user
TARGET_USER="${1:-azureuser}"
echo "source /usr/local/gromacs/bin/GMXRC" >> /home/$TARGET_USER/.bashrc
chown $TARGET_USER:$TARGET_USER /home/$TARGET_USER/.bashrc
