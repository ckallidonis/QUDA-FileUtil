#!/bin/bash

if [ $# -ne 4 ]; 
    then echo "Usage: $0 <inifile_base> <outfile_base> <mom_list> <N_GPU>"
    exit
fi

INI_BASE=$1
OUT_BASE=$2
MOM_LIST=$3
NGPU=$4

for tp in Loops LoopsCv LpsDw LpsDwCv
do
    OUT_FILE=${OUT_BASE}_rshp_${tp}.loop.dat
    rm -f ${OUT_FILE}

    while read x y z
    do
	for dir in 00 01 02 03
	do
	    for((i=0;i<${NGPU};i++)){
		FILE=${INI_BASE}_${tp}.loop.${NGPU}_${i}
		
		grep -- "^.. .. ${dir} ${x} ${y} ${z}" ${FILE} >> ${OUT_FILE}
	    }
	done
    done < ${MOM_LIST}
    echo $tp
done


for tp in Scalar dOp
do
    OUT_FILE=${OUT_BASE}_rshp_${tp}.loop.dat
    rm -f ${OUT_FILE}

    while read x y z
    do
	for((i=0;i<${NGPU};i++)){
	    FILE=${INI_BASE}_${tp}.loop.${NGPU}_${i}
	    
	    grep -- "^.. .. ${x} ${y} ${z}" ${FILE} >> ${OUT_FILE}
	}
    done < ${MOM_LIST}

    echo $tp
done
