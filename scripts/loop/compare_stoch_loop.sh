#!/bin/bash

if [ $# -ne 7 ]; 
    then echo "Usage: $0 <loop_base_1> <loop_base_2> <stoch_list> <paste_dir> <compare_dir> <type 0:Local, 1:oneD> <remove-matrix-elements 0,1>"
    exit
fi

LOOP_BASE1=$1
LOOP_BASE2=$2
NSTOCH_LIST=$3
PASTE_BASE=$4/paste
COMPR_BASE=$5/compare
tp=$6
melem=$7

if [ ${tp} -lt 0 ]
then
    echo "type must be either 0:Local or 1:oneD"
    exit
fi
if [ ${tp} -gt 1 ]
then
    echo "type must be either 0:Local or 1:oneD"
    exit
fi

if [ ${tp} -eq 0 ]
then

    for ctype in Scalar dOp
    do
	for Ns in `cat $NSTOCH_LIST`
	do
	    LOOP_FILE1=${LOOP_BASE1}_${ctype}.loop.${Ns}.dat
	    LOOP_FILE2=${LOOP_BASE2}_${ctype}.loop.${Ns}.dat
	    PASTE_FILE=${PASTE_BASE}_${ctype}.loop.${Ns}.dat
	    COMPR_FILE=${COMPR_BASE}_${ctype}.loop.${Ns}.dat

	    rm -f $PASTE_FILE $COMPR_FILE

	    paste ${LOOP_FILE1} ${LOOP_FILE2} > ${PASTE_FILE}
	    
	    cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %8.6f %8.6f\n",$1,$2,$6/$13,$7/$14)}' > ${COMPR_FILE}

	    if [ ${melem} -eq 1 ]
	    then
		COMPR_FILE1=${COMPR_BASE}_${ctype}.loop.${Ns}.temp
		cat ${COMPR_FILE} | grep -v "^.. 02" | grep -v "^.. 07" | grep -v "^.. 08" | grep -v "^.. 13" > ${COMPR_FILE1}		
		mv ${COMPR_FILE1} ${COMPR_FILE}
	    fi

	    nl=`wc -l ${COMPR_FILE}  | awk '{print $1}'`
		    
	    echo "$ctype , $Ns: N = $nl"

	    cat ${COMPR_FILE} | awk '{s1+=$3} {s2+=$4} END {printf("                                %.4f   %.4f\n",s1,s2)}'
	done
    done
fi

if [ ${tp} -eq 1 ]
then

    for ctype in Loops LoopsCv LpsDw LpsDwCv
    do
	for Ns in `cat $NSTOCH_LIST`
	do
	    LOOP_FILE1=${LOOP_BASE1}_${ctype}.loop.${Ns}.dat
	    LOOP_FILE2=${LOOP_BASE2}_${ctype}.loop.${Ns}.dat
	    PASTE_FILE=${PASTE_BASE}_${ctype}.loop.${Ns}.dat
	    COMPR_FILE=${COMPR_BASE}_${ctype}.loop.${Ns}.dat
	    
	    rm -f $PASTE_FILE $COMPR_FILE
	    
	    paste ${LOOP_FILE1} ${LOOP_FILE2} > ${PASTE_FILE}
	    
	    cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %02.0f %8.6f %8.6f\n",$1,$2,$3,$7/$15,$8/$16)}' > ${COMPR_FILE}
	    
	    if [ ${melem} -eq 1 ]
	    then
		COMPR_FILE1=${COMPR_BASE}_${ctype}.loop.${Ns}.temp
		cat ${COMPR_FILE} | grep -v "^.. 02" | grep -v "^.. 07" | grep -v "^.. 08" | grep -v "^.. 13" > ${COMPR_FILE1}		
		mv ${COMPR_FILE1} ${COMPR_FILE}
	    fi

	    nl=`wc -l ${COMPR_FILE}  | awk '{print $1}'`
		    
	    echo "$ctype , $Ns: N = $nl"

	    cat ${COMPR_FILE} | awk '{s1+=$4} {s2+=$5} END {printf("                                %.4f   %.4f\n",s1,s2)}'
	done
    done
fi
