#!/bin/sh
debug1: client_input_channel_req: channel 0 rtype keepalive@openssh.com reply 1
module unload  Intel/12.0
module unload  Intel/13.0
module unload Intel
module load Intel/15.x
module load Openmpi/intel-1.8.4
module load CMake/3.2.2
module load HDF5/1.8.8
HDF5_ROOT=/export/zgid/zgid/software/hdf5-1.8.10/lib
module unload  Intel/13.x
module unload  Intel/13.x
CC=mpicc FC=mpif90 HDF5_ROOT="/export/zgid/zgid/software/hdf5-1.8.10/" cmake ../
make -j 8


#ifort -O3 -g -openmp *.f90 -mkl=sequential
