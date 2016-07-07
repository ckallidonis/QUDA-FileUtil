#!/bin/bash

if [ $# -ne 8 ]; 
    then echo "Usage: $0 <src_list> <threep_base_1> <threep_base_2> <paste_dir> <compare_dir> <tsink> <proj> <type 0:noether, 1:ultra_local, 2:oneD>"
    exit
fi

SRC_LIST=${1}
THRP_BASE1=${2}
THRP_BASE2=${3}
PASTE_BASE=${4}/paste
COMPR_BASE=${5}/compare
TSINK=$6
proj=$7
tp=$8

if [ ${tp} -lt 0 ]
then
    echo "type must be either 0:noether, 1:ultra_local, 2:oneD"
    exit
fi
if [ ${tp} -gt 2 ]
then
    echo "type must be either 0:noether, 1:ultra_local, 2:oneD"
    exit
fi
if [ ${tp} -eq 0 ]
then
    ctype=noether
fi
if [ ${tp} -eq 1 ]
then
    ctype=ultra_local
fi
if [ ${tp} -eq 2 ]
then
    ctype=oneD
fi

echo "Running for type ${ctype}"

while read x y z t
do
    src=${x}.${y}.${z}.${t}

    for tsink in $TSINK
    do
	for part in up down
	do
	    THRP_FILE1=${THRP_BASE1}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${src}.dat
	    THRP_FILE2=${THRP_BASE2}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${src}.dat
	    PASTE_FILE=${PASTE_BASE}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${src}.dat
	    COMPR_FILE=${COMPR_BASE}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${src}.dat
	    
	    rm -f ${PASTE_FILE} ${COMPR_FILE}

	    paste ${THRP_FILE1} ${THRP_FILE2} > ${PASTE_FILE}
	    
	    if [ ${tp} -eq 2 ]
	    then
		cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %8.5f %8.5f\n",$1,$2,$7/$15,$8/$16)}' > ${COMPR_FILE}
	    fi
	    
	    if [ ${tp} -lt 2 ]
	    then
		cat ${PASTE_FILE}  | awk '{printf("%02.0f %02.0f %8.5f %8.5f\n",$1,$2,$6/$13,$7/$14)}' > ${COMPR_FILE}
	    fi

	    nl=`wc -l ${COMPR_FILE}  | awk '{print $1}'`
	    
	    echo "$src , $tsink , $part: N = $nl"
	    cat ${COMPR_FILE} | awk '{s1+=$3} {s2+=$4} END {printf("                                %.4f   %.4f\n",s1,s2)}'

	done
    done    
done < ${SRC_LIST}
