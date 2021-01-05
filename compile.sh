#!/bin/bash
set -x
# B. Pabla Oct/2019
# script to compile Joana's code using netcdf for i/o 
#loads intel compiler and MPI library ...mandatory for eccc-ppp3,4
. ssmuse-sh -x comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199
. ssmuse-sh -x hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.4.0-mofed-4.6--intel-19.0.3.199
. ssmuse-sh -x hpco/exp/openmpi-setup/openmpi-setup-0.1

# load s.compile
. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/01.3
# load rmnlib
. r.load.dot rpn/libs/19.4
. r.load.dot rpn/utils/19.4
# loads netcdf (serial)stuff
#. ssmuse-sh -x hpco/exp/hdf5-netcdf4/serial/static/intel-19.0.3.199/02
# For parallel version of HDF5/NetCDF
. ssmuse-sh -x hpco/exp/hdf5-netcdf4/parallel/openmpi-3.1.2/static/intel-19.0.3.199/01
#
#s.compile -src kz_hc_GEM-MACH_v24_read_clust_v3.f90 -o joana.abs -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc
s.compile -src kz_hc_GEM-MACH_v24_read_clust_v3.f90 -o kzCluster.abs -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc
#s.compile -src kz_hc_GEM-MACH_v24_read_clust_v3.f90 -o joana.abs -optf="-C -pg -traceback -fp-model strict -mp1 -warn all" -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc
