#!/bin/bash
set -x

. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2

. ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2022.1.2/shared/inteloneapi-2022.1.2/01

. r.load.dot rpn/libs/20220216
. r.load.dot rpn/utils/20220408
. r.load.dot rpn/vgrid/20220216
. r.load.dot rpn/ezinterpv/20220216

# s.compile -src mheap.f90 -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -O 1 # -debug -strict -prof -O 0
s.compile -src mheap.f90 -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -debug -strict -prof -O 0

#s.compile -src cluster_mpi.f90 -openmp -mpi -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -O 1 
s.compile -src cluster_mpi.f90 -openmp -mpi -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -debug -strict -prof -O 0

# s.compile -obj cluster_mpi.o mheap.o -o cluster.abs -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -O 1 # -debug -strict -prof -O 0
s.compile -obj cluster_mpi.o mheap.o -o cluster.abs -openmp -mpi  -librmn -libappl netcdff netcdf hdf5_hl hdf5 z curl irc -debug -strict -prof -O 0
