#!/bin/bash

if [ $# -ne 5 ]
then
    echo "Usage: $0 <extract_dir> <out_dir> <mom_list> <stoch_list> <traj>"
    exit
fi

EDIR=$1
OUT_DIR=$2
MOM_LIST=$3
NSTOCH_LIST=$4
TRAJ=$5

dir=${EDIR}/stoch_hdf5_extracted

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    for Ns in `cat $NSTOCH_LIST`
    do
	OUT_FILE=${OUT_DIR}/loop_stoch.${TRAJ}_${tp}.loop.${Ns}.dat

	while read mx my mz
	do
	    FILE=${dir}/loop_stoch.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.${Ns}.dat

	    cat ${FILE}

	done < ${MOM_LIST} > ${OUT_FILE}
    done

    echo "${tp}"
done
