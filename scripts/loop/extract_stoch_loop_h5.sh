#!/bin/bash

if [ $# -ne 8 ]
then
    echo "Usage: $0 <h5-file> <extract_dir> <out_dir> <mom_list> <stoch_list> <exe_dir> <T> <traj>"
    exit
fi


FILE=$1
DIR=$2
OUT_DIR=$3
MOM_LIST=$4
NSTOCH_LIST=$5
EXE_DIR=$6
T=$7
TRAJ=$8

EX_DIR=${DIR}/stoch_hdf5_extracted

mkdir -p ${EX_DIR}

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    for Ns in `cat $NSTOCH_LIST`
    do
	OUT_FILE=${OUT_DIR}/loop_stoch.${TRAJ}_${tp}.loop.${Ns}.dat
	rm -f ${OUT_FILE}

	while read mx my mz
	do
	    EX_FILE=${EX_DIR}/loop_stoch.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.${Ns}.dat
	    
	    ${EXE_DIR}/extract_stoch_loop_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${TRAJ} ${Ns} ${T}

	    cat ${EX_FILE} >> ${OUT_FILE}
	done < ${MOM_LIST}
    done
done
