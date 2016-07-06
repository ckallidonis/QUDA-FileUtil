#!/bin/bash

if [ $# -ne 9 ]
then
    echo "Usage: $0 <ini_dir> <out_dir> <src_list> <mom_list> <conf> <tsink> <proj> <T> <type 0:noether, 1:ultra_local, 2:oneD>"
    exit
fi

INI_DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
TSINK=$6
PROJ=$7
T=$8
tp=$9

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


while read x y z t
do
    for part in up down
    do
	for tsink in $TSINK
	do
	    for proj in $PROJ
	    do

		THRP_FILE=${INI_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${x}.${y}.${z}.${t}.dat
		OUT_FILE=${OUT_DIR}/threep_rshp.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${ctype}.SS.${x}.${y}.${z}.${t}.dat

		while read mx my mz
		do
		    if [ ${tp} -ne 2 ]
                    then
			for((ts=0;ts<=${tsink};ts++)){
			    grep -- "${mx} ${my} ${mz}" ${THRP_FILE} | awk '(NR-1)%'${T}' == '${ts}''
			} 
		    fi
		    
		    if [ ${tp} -eq 2 ]
                    then
			for((ts=0;ts<=${tsink};ts++)){ grep -- "${mx} ${my} ${mz}" ${THRP_FILE} | awk '(NR-1)%'${T}' == '${ts}''; } | sort -n -k2 -k3 -k1
		    fi
		done < $MOM_LIST > ${OUT_FILE}
	    done
       	done
    done
    
    echo sx${x}sy${y}sz${z}st${t}

done < $SRC_LIST
