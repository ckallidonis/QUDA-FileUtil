#!/bin/bash

if [ $# -ne 7 ]
then
    echo "Usage: $0 <extract_dir> <out_dir> <src_list> <mom_list> <conf> <tsink> <proj>"
    exit
fi

DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
TSINK=$6
PROJ=$7

while read x y z t
do
    src=sx${x}sy${y}sz${z}st${t}
    EX_DIR=${DIR}/${src}

    for tsink in $TSINK
    do
	for proj in $PROJ
	do
	    for part in up down
	    do
		for tp in ultra_local noether oneD
		do
		    OUT_FILE=${OUT_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.SS.${x}.${y}.${z}.${t}.dat
		    
		    while read mx my mz
		    do
			FILE=${EX_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
			
			cat ${FILE}
			
		    done < ${MOM_LIST} > ${OUT_FILE}
		done
	    done
	done
    done
    
    echo $src

done < ${SRC_LIST}
