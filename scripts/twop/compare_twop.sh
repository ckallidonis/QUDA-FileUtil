#!/bin/bash

if [ $# -ne 6 ]; 
    then echo "Usage: $0 <src_list> <twop_base_1> <twop_base_2> <paste_dir> <compare_dir> <type 0:mesons, 1:baryons>"
    exit
fi

SRC_LIST=${1}
TWOP_BASE1=${2}
TWOP_BASE2=${3}
PASTE_BASE=${4}/paste
COMPR_BASE=${5}/compare
tp=$6
Nth=33


if [ ${tp} -lt 0 ]
then
    echo "type must be either 0:mesons or 1:baryons"
    exit
fi
if [ ${tp} -gt 1 ]
then
    echo "type must be either 0:mesons or 1:baryons"
    exit
fi
if [ ${tp} -eq 0 ]
then
    ctype=mesons
fi
if [ ${tp} -eq 1 ]
then
    ctype=baryons
fi


while read x y z t
do
    src=${x}.${y}.${z}.${t}

    TWOP_FILE1=${TWOP_BASE1}.${ctype}.SS.${src}.dat
    TWOP_FILE2=${TWOP_BASE2}.${ctype}.SS.${src}.dat
    PASTE_FILE=${PASTE_BASE}.${ctype}.SS.${src}.dat
    COMPR_FILE=${COMPR_BASE}.${ctype}.SS.${src}.dat

    rm -f ${PASTE_FILE} ${COMPR_FILE}


    paste ${TWOP_FILE1} ${TWOP_FILE2} > ${PASTE_FILE}

    if [ ${tp} -eq 1 ]
    then
	cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %8.6f %8.6f %8.6f %8.6f\n",$1,$2,$8/$19,$9/$20,$10/$21,$11/$22)}' > ${COMPR_FILE}

    fi

    if [ ${tp} -eq 0 ]
    then
	COMPR_FILE1=${COMPR_BASE}_all.${ctype}.SS.${src}.dat
	COMPR_FILE2=${COMPR_BASE}_No1stline.${ctype}.SS.${src}.dat

	cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %8.6f %8.6f %8.6f %8.6f\n",$1,$2,$6/$15,$7/$16,$8/$17,$9/$18)}' > ${COMPR_FILE1}

	tail -n +2 ${COMPR_FILE1} > ${COMPR_FILE2}

	awk 'NR % 33 != 0' ${COMPR_FILE2} > ${COMPR_FILE}
    fi

    nl=`wc -l ${COMPR_FILE}  | awk '{print $1}'`

    echo "$src: N = ${nl}"
    cat ${COMPR_FILE} | awk '{s1+=$3} {s2+=$4} {s3+=$5} {s4+=$6} END {printf("%.4f   %.4f    %.4f   %.4f\n",s1,s2,s3,s4)}'
        
done < ${SRC_LIST}
