#!/bin/bash

if [ $# -ne 7 ]
then
    echo "Usage: $0 <h5-file> <extract_dir> <out_dir> <mom_list> <exe_dir> <T> <traj>"
    exit
fi


FILE=$1
DIR=$2
OUT_DIR=$3
MOM_LIST=$4
EXE_DIR=$5
T=$6
TRAJ=$7

EX_DIR=${DIR}/exact_hdf5_extracted

mkdir -p ${EX_DIR}

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    OUT_FILE=${OUT_DIR}/loop_exact.${TRAJ}_${tp}.loop.dat
    rm -f ${OUT_FILE}

    while read mx my mz
    do
	EX_FILE=${EX_DIR}/loop_exact.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.dat
	    
	${EXE_DIR}/extract_exact_loop_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${TRAJ} ${T}

	cat ${EX_FILE} >> ${OUT_FILE}
    done < ${MOM_LIST}
done

rm -rf ${EX_DIR}
