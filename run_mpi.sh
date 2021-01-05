#!/bin/bash
rm joana_batch.scr
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

r.run_in_parallel -npex 40 -inorder -pgm \${wrkDir}/cluster.abs
EOD
# on ppp3 -cpus 40x1 -cm 5G gives only one node; as a workaround use -cpus 44x44 -cm 210G
ord_soumet jobscript.sh -mach eccc-ppp3 -cpus 40x40 -cm 100G -waste 90 -w 15 -mpi -jn colin_parallel -listing /space/hall3/sitestore/eccc/aq/r1/cle001/src/parallel/listings
#ord_soumet joana_batch.scr -mach eccc-ppp4 -cpus 4x1 -cm 10G -waste 90 -t 16200 -mpi -jn pabla -listing /home/bap001/joana_networkanalysis/v5_netcdf/parallel/listings

