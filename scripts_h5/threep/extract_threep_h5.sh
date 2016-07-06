#!/bin/bash

if [ $# -ne 8 ]
then
    echo "Usage: $0 <extact_dir> <src_list> <mom_list> <conf> <tsink> <proj> <Qsq> <exe_dir>"
    exit
fi

DIR=$1
SRC_LIST=$2
MOM_LIST=$3
CONF=$4
TSINK=$5
PROJ=$6
Qsq=$7
EXE_DIR=$8

while read x y z t
do 
    FILE=${DIR}/threep.${CONF}_neutron_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}

    OUT_DIR=${DIR}/${src}

    mkdir -p ${OUT_DIR}

    for tsink in $TSINK
    do 
	for proj in $PROJ
	do
	    for part in up down
	    do
		for tp in ultra_local noether oneD
		do
		    while read mx my mz
		    do
			OUT_FILE=${OUT_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
			
			${EXE_DIR}/extract_threep_h5 ${FILE} ${OUT_FILE} ${tp} ${proj} ${part} ${mx} ${my} ${mz} ${CONF} ${src} ${tsink}
		    done < ${MOM_LIST}
		done
	    done
	done
    done
done < ${SRC_LIST}
