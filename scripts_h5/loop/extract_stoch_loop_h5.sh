#!/bin/bash

if [ $# -ne 7 ]
then
    echo "Usage: $0 <h5-file> <extract_dir> <mom_list> <stoch_list> <exe_dir> <T> <traj>"
    exit
fi


FILE=$1
EDIR=$2
MOM_LIST=$3
NSTOCH_LIST=$4
EXE_DIR=$5
T=$6
TRAJ=$7

dir=${EDIR}/stoch_hdf5_extracted

mkdir -p ${dir}

for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
do
    for Ns in `cat $NSTOCH_LIST`
    do
	while read mx my mz
	do
	    OUT_FILE=${dir}/loop_stoch.${TRAJ}_${tp}.${mx}_${my}_${mz}.loop.${Ns}.dat
	    
	    ${EXE_DIR}/extract_stoch_loop_h5 ${FILE} ${OUT_FILE} ${tp} ${mx} ${my} ${mz} ${TRAJ} ${Ns} ${T}
	done < ${MOM_LIST}
    done
done
