#!/bin/bash

# set up batch system parameters for MPI threads and wall clock time (min)
# this setup should be suitable for a useFract = 6 or up
NODESX=16   # number of nodes
NODESY=1   # MPI threads per node
WTIME=270    # wall clock time in minutes
# use this setup for using all the data useFract=1
# NODESX=40   # number of nodes
# NODESY=1   # MPI threads per node
# WTIME=60    # wall clock time in minutes

# Assume we want to use a whole node so set OMPTHREADS and MEM
# to the cores and memory per node divided by the number of MPI
# threads per node. OMPTHREADS is not the same as the OMP variable
# OMP_NUM_THREADS because that should be set by the grid engine.
OMPTHREADS=$(( 40 / NODESY ))
MEM="$(( 160 / NODESY ))G"

rm jobscript.sh
cat <<EOD >jobscript.sh
#loads intel compiler and MPI library ...mandatory for eccc-ppp3,4
. ssmuse-sh -x comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199
. ssmuse-sh -x hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.4.0-mofed-4.6--intel-19.0.3.199
. ssmuse-sh -x hpco/exp/openmpi-setup/openmpi-setup-0.1

# load s.compile
. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/01.3
# load rmnlib
. r.load.dot rpn/libs/19.4
. r.load.dot rpn/utils/19.4

# For parallel version of HDF5/NetCDF
. ssmuse-sh -x hpco/exp/hdf5-netcdf4/parallel/openmpi-3.1.2/static/intel-19.0.3.199/01
#

wrkDir=/space/hall3/sitestore/eccc/aq/r1/cle001/src/parallel
cd \${wrkDir}

# export OMPI_MCA_io_base_verbose=40

r.run_in_parallel -npex ${NODESX} -npey ${NODESY} -inorder -pgm \${wrkDir}/cluster.abs -verbose -args \$@
EOD
# on ppp3 -cpus 40x1 -cm 5G gives only one node; as a workaround use -cpus 44x44 -cm 210G

ord_soumet jobscript.sh -mach eccc-ppp3 -cpus ${NODESX}x${NODESY}x${OMPTHREADS} -cm ${MEM} -waste 90 -w ${WTIME} -mpi -jn colin_mpi -listing /space/hall3/sitestore/eccc/aq/r1/cle001/src/parallel/listings -args $@
