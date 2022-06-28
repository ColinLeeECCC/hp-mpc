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
    # Get the MPI configuration
    node_list=$( sed -n -e "s/^.*Initial node list = \(.*\)$/\1/p" $f | tr ' ' '\n' | uniq -c)
    mpi_x=$( echo "$node_list" | wc -l )
    mpi_y=$( echo "$node_list" | awk '{ print $1 }' | uniq )
    dims=$( sed -n -e "s/^oe-00000:  dims *\([0-9]*\) *\([0-9]*\) *1.*$/\1,\2/p" $f )
    d1=${dims%,*}
    d2=${dims#*,}
    size=$( expr $d1 \* $d2 )
    omp_threads=$( sed -n -e "s/^.*++ OMP_NUM_THREADS=\([0-9]\+\)/\1/p" $f )
    #size=$d1
    numdays=$( grep "^oe-00000:  Loading data for day" $f | wc -l )
    ## I'm not sure why this is true, but the oe-##### seems to be linked to the node maybe? So there are several
    ## oe-00000: groups in the code, so we need to divide numdays by mpi_y
    ## Maybe this only happens on the frontend machines...
    if [[ $f == *"ppp5-"* || $f == *"ppp6-"* ]]; then
	echo "$f -> dividing numdays"
	numdays=$(( numdays / mpi_y ))
    fi
    
    t1=$( sed -n -e "s/^oe-00000:   Calcu\?lating dissimilarity matrix took *\([0-9E\.\-]*\) *seconds/\1/p" $f | head -n 1 )
    t2=$( sed -n -e "s/^oe-00000:  Loading dissimilarity matrix took *\([0-9E\.\-]*\) *seconds/\1/p" $f | head -n 1)
    t3=$( sed -n -e "s/^oe-00000:  Clust\?ering took *\([0-9E\.\-]*\) *sec/\1/p" $f | head -n 1 ) # turns out there's a typo in the code

    ttot=$( sed -n -e "s/^.*WallTime=\([0-9a-z]\+\) .*$/\1/p" $f )

    if [[ "$ttot" == *"s" ]]; then
	ttot=${ttot%s}
    elif [[ "$ttot" == *"m" ]]; then
	ttot=${ttot%m}
	tott=$((ttot*60))
    elif [[ "$ttot" == *"h" ]]; then
	ttot=${ttot%h}
	tott=$((ttot*60))
    else
	[[ -z "$ttot" ]] || echo "Couldn't figure out time units for total WallTime"
    fi
    

    if [[ -z "$mpi_y" ]]; then
	echo "Something is wrong with $f" >&2
	shift
	continue
    fi
    # numdays=$(( numdays / mpi_y ))
    
    [ -z "$t1" ] && t1="0"
    [ -z "$t2" ] && t2="0"
    [ -z "$t3" ] && t3="0"

    echo -e "$size\t$mpi_x\t$mpi_y\t$omp_threads\t$numdays\t$t1\t$t2\t$t3\t$ttot"
    shift
done
    
