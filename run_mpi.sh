#!/bin/bash

# set up batch system parameters for MPI threads and wall clock time (min)
# this setup should be suitable for a useFract = 6 or up
NODESX=10     # number of nodes
NODESY=8    # MPI threads per node
WTIME=180    # wall clock time in minutes
# use this setup for using all the data useFract=1
# NODESX=40   # number of nodes
# NODESY=1   # MPI threads per node
# WTIME=60    # wall clock time in minutes

# Assume we want to use a whole node so set OMPTHREADS and MEM
# to the cores and memory per node divided by the number of MPI
# threads per node. OMPTHREADS is not the same as the OMP variable
# OMP_NUM_THREADS because that should be set by the grid engine.
# OMPTHREADS=$(( 80 / NODESY ))
# MEM="$(( 100 / NODESY ))G"
OMPTHREADS=10
MEM=16G

rm jobscript.sh
cat <<EOD >jobscript.sh
. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2

. ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2022.1.2/shared/inteloneapi-2022.1.2/01

. r.load.dot rpn/libs/20220216
. r.load.dot rpn/utils/20220408
. r.load.dot rpn/vgrid/20220216
. r.load.dot rpn/ezinterpv/20220216

wrkDir=/space/hall5/sitestore/eccc/aq/r1/cle001/src/hierarchical-clustering
cd \${wrkDir}

# export OMPI_MCA_io_base_verbose=40

r.run_in_parallel -npex ${NODESX} -npey ${NODESY} -inorder -pgm \${wrkDir}/cluster.abs -verbose -args \$@
EOD
# on ppp3 -cpus 40x1 -cm 5G gives only one node; as a workaround use -cpus 44x44 -cm 210G

ord_soumet jobscript.sh -mach ppp6 -cpus ${NODESX}x${NODESY}x${OMPTHREADS} -cm ${MEM} -waste 90 -w ${WTIME} -mpi -jn mpi_test -listing /space/hall5/sitestore/eccc/aq/r1/cle001/src/hierarchical-clustering/listings -args $@
