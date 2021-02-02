#!/bin/bash

# extract timing info from listings from cluster_omp.abs runs
set +x

USAGE="Usage: $0 <listings>"

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

while (( "$#" )); do
    f=$1
    # echo $f
    dims=$( sed -n -e "s/^ dims *\([0-9]*\) *\([0-9]*\) *1.*$/\1,\2/p" $f )
    d1=${dims%,*}
    d2=${dims#*,}
    size=$( expr $d1 \* $d2 )
    #size=$d1
    numdays=$( grep "^ Loading data for day" $f | wc -l )
    t1=$( sed -n -e "s/^ Loading data and precalculating took *\([0-9E\.\-]*\) *sec/\1/p" $f )
    t2=$( sed -n -e "s/^ Calculating initial dissimilarity scores took *\([0-9E\.\-]*\) *sec/\1/p" $f )
    t3=$( sed -n -e "s/^ Clust\?ering took *\([0-9E\.\-]*\) *sec/\1/p" $f ) # turns out there's a typo in the code
    [ -z "$t1" ] && t1="0"
    [ -z "$t2" ] && t2="0"
    [ -z "$t3" ] && t3="0"

    echo "$size $numdays $t1 $t2 $t3"
    shift
done
    
