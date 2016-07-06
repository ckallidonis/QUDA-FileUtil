#!/bin/bash

if [ $# -ne 6 ]
then
    echo "Usage: $0 <h5-file> <extract_dir> <mom_list> <exe_dir> <T> <traj>"
    exit
fi


FILE=$1
EDIR=$2
MOM_LIST=$3
EXE_DIR=$4
T=$5
TRAJ=$6

dir=${EDIR}/exact_hdf5_extracted

mkdir -p ${dir}

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    while read mx my mz
    do
	OUT_FILE=${dir}/loop_exact.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.dat
	    
	${EXE_DIR}/extract_exact_loop_h5 ${FILE} ${OUT_FILE} ${tp} ${mx} ${my} ${mz} ${TRAJ} ${T}
    done < ${MOM_LIST}
done
