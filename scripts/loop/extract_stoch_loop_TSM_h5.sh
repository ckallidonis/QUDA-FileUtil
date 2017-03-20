#!/bin/bash

if [ $# -ne 11 ]
then
    echo "Usage: $0 <h5-file> <extract_dir> <out_dir> <mom_list> <stoch_list> <LowPrecSum:1, HighPrecSum:0> <HighPrec, LowPrec> <exe_dir> <T> <traj> <HighMomForm:1 LowMomForm:0>"
    exit
fi


FILE=$1
DIR=$2
OUT_DIR=$3
MOM_LIST=$4
NSTOCH_LIST=$5
SUM_TYPE=$6
PREC_TYPE=$7
EXE_DIR=$8
T=$9
TRAJ=${10}
HighMomForm=${11}

EX_DIR=${DIR}/stoch_hdf5_extracted

mkdir -p ${OUT_DIR}

if [ $HighMomForm -eq 0 ]
then

    mkdir -p ${EX_DIR}

    if [ ${SUM_TYPE} -eq 1 ]
    then
	STRING=NLP
    fi
    if [ ${SUM_TYPE} -eq 0 ]
    then
	STRING=${PREC_TYPE}_NHP
    fi

    for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
    do
	for Ns in `cat $NSTOCH_LIST`
	do
            OUT_FILE=${OUT_DIR}/loop_stoch.${TRAJ}_TSM_${STRING}${Ns}_${tp}.loop.dat
	    rm -f ${OUT_FILE}

	    while read mx my mz
	    do
		EX_FILE=${EX_DIR}/loop_stoch.${TRAJ}_${STRING}${Ns}_${tp}.${mx}_${my}_${mz}.loop.dat
	        
		${EXE_DIR}/extract_stoch_loop_TSM_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${TRAJ} ${Ns} ${SUM_TYPE} ${T}

		cat ${EX_FILE} >> ${OUT_FILE}
	    done < ${MOM_LIST}
	done
    done

    rm -rf ${EX_DIR}

fi


if [ $HighMomForm -eq 1 ]
then

    if [ ${SUM_TYPE} -eq 1 ]
    then
	STRING=NLP
    fi
    if [ ${SUM_TYPE} -eq 0 ]
    then
	STRING=${PREC_TYPE}_NHP
    fi

    for tp in Scalar dOp Loops LoopsCv LpsDw LpsDwCv
    do
	for Ns in `cat $NSTOCH_LIST`
	do
            OUT_FILE=${OUT_DIR}/loop_stoch.${TRAJ}_TSM_${STRING}${Ns}_${tp}.loop.dat
	    rm -f ${OUT_FILE}
	        
	    ${EXE_DIR}/extract_stoch_loop_TSM_HighMomForm_h5 ${FILE} ${OUT_FILE} ${tp} ${TRAJ} ${Ns} ${SUM_TYPE} ${T}
		
	done
    done

fi
