#!/bin/bash

if [ $# -ne 5 ]
then
    echo "Usage: $0 <ini_dir> <out_dir> <src_list> <mom_list> <conf>"
    exit
fi

INI_DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5

while read x y z t
do
    MES_FILE=${INI_DIR}/twop.${CONF}.mesons.SS.${x}.${y}.${z}.${t}.dat
    MESOUT_FILE=${OUT_DIR}/twop_rshp.${CONF}.mesons.SS.${x}.${y}.${z}.${t}.dat
    
    rm -f ${MESOUT_FILE}

    for i in 0 1 2 3 4 5 6 7 8 9
    do
	while read mx my mz
	do
	    grep -- "^${i}" ${MES_FILE} | grep -- "${mx} ${my} ${mz}" >> ${MESOUT_FILE}
	done < ${MOM_LIST}
    done

    echo sx${x}sy${y}sz${z}st${t}

done < ${SRC_LIST}
