#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Usage: $0 <extract_dir> <out_dir> <mom_list> <traj>"
    exit
fi

EDIR=$1
OUT_DIR=$2
MOM_LIST=$3
TRAJ=$4

dir=${EDIR}/exact_hdf5_extracted

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    OUT_FILE=${OUT_DIR}/loop_exact.${TRAJ}_${tp}.loop.dat

    while read mx my mz
    do
	FILE=${dir}/loop_exact.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.dat
	
	cat ${FILE}
	
    done < ${MOM_LIST} > ${OUT_FILE}

    echo "${tp}"
done
