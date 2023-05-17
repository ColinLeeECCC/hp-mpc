#!/bin/bash

# set -x

eval $(cclargs_lite -D , $0     "[test cclargs]"\
 -nodesx "1" ""   "[number of nodes to use]"\
 -nodesy "8" ""   "[number of cores per node]"\
 -ompnum "10" ""  "[number of openMP threads]"\
 -wtime  "270" "" "[wallclock time in minutes]"\
 -mem    "" ""  "[memory requested]"\
 ++ $*)

NODESX="$nodesx"
NODESY="$nodesy"
OMPTHREADS="$ompnum"
WTIME="$wtime"
defaultmem=$( bc <<< "scale=0; 160000 * ${ompnum} / ${nodesy} / 80 " )
MEM="${mem:-${defaultmem}M}"

rm jobscript.sh
cat <<EOD >jobscript.sh
. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2

. ssmuse-sh -d main/opt/hdf5-netcdf4/parallel/intelmpi-2022.1.2/shared/inteloneapi-2022.1.2/01

. r.load.dot rpn/libs/20220216
. r.load.dot rpn/utils/20220408
. r.load.dot rpn/vgrid/20220216
. r.load.dot rpn/ezinterpv/20220216

wrkDir=/space/hall5/sitestore/eccc/aq/r1/cle001/src/hierarchical-clustering2
cd \${wrkDir}

# export OMPI_MCA_io_base_verbose=40

r.run_in_parallel -npex ${NODESX} -npey ${NODESY} -inorder -pgm \${wrkDir}/cluster.abs -verbose -args \$@
EOD
# on ppp3 -cpus 40x1 -cm 5G gives only one node; as a workaround use -cpus 44x44 -cm 210G

echo "NODESX=$NODESX  NODESY=$NODESY  OMPTHREADS=$OMPTHREADS"
echo "WTIME=$WTIME  MEM=$MEM  ARGS=$@"
ord_soumet jobscript.sh -mach ppp5 -cpus ${NODESX}x${NODESY}x${OMPTHREADS} -cm ${MEM} -waste 99 -w ${WTIME} -mpi -jn hc_bulk_timing -listing /space/hall5/sitestore/eccc/aq/r1/cle001/src/hierarchical-clustering2/listings -args $@
