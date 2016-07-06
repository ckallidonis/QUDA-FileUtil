#!/bin/bash

if [ $# -ne 5 ]; 
    then echo "Usage: $0 <inifile_base> <outfile_base> <mom_list> <stoch_list> <N_GPU>"
    exit
fi

INI_BASE=$1
OUT_BASE=$2
MOM_LIST=$3
NSTOCH_LIST=$4
NGPU=$5

for tp in Loops LoopsCv LpsDw LpsDwCv
do
    for Ns in `cat ${NSTOCH_LIST}`
    do
	OUT_FILE=${OUT_BASE}_rshp_${tp}.loop.${Ns}.dat
	rm -f ${OUT_FILE}

	while read x y z
	do
	    for dir in 00 01 02 03
	    do
		for((i=0;i<${NGPU};i++)){
		    FILE=${INI_BASE}_${tp}.loop.${Ns}.${NGPU}_${i}
		    
		    grep -- "^.. .. ${dir} ${x} ${y} ${z}" ${FILE} >> ${OUT_FILE}
		}
	    done
	done < $MOM_LIST
    done

    echo $tp
done


for tp in Scalar dOp
do
    for Ns in `cat ${NSTOCH_LIST}`
    do
	OUT_FILE=${OUT_BASE}_rshp_${tp}.loop.${Ns}.dat
	rm -f ${OUT_FILE}

	while read x y z
	do
	    for((i=0;i<${NGPU};i++)){
		FILE=${INI_BASE}_${tp}.loop.${Ns}.${NGPU}_${i}
		
		grep -- "^.. .. ${x} ${y} ${z}" ${FILE} >> ${OUT_FILE}
	    }
	done < $MOM_LIST
    done

    echo $tp
done
