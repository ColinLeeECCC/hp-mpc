#!/bin/bash

. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2

. ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2022.1.2/shared/inteloneapi-2022.1.2/01

. r.load.dot rpn/libs/20220216

s.compile -src csv_to_netcdf.F90 -o csv_to_netcdf.abs -mpi -libappl netcdff netcdf hdf5_hl hdf5 z curl irc
