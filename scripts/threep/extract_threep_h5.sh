#!/bin/bash

if [ $# -ne 9 ]
then
    echo "Usage: $0 <extact_dir> <out_dir> <src_list> <mom_list> <conf> <tsink> <proj> <Qsq> <exe_dir>"
    exit
fi

DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
TSINK=$6
PROJ=$7
Qsq=$8
EXE_DIR=$9

while read x y z t
do 
    FILE=${DIR}/threep.${CONF}_neutron_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}

    EX_DIR=${DIR}/${src}

    mkdir -p ${EX_DIR}

    for tsink in $TSINK
    do 
	for proj in $PROJ
	do
	    for part in up down
	    do
		for tp in ultra_local noether oneD
		do
		    OUT_FILE=${OUT_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.SS.${x}.${y}.${z}.${t}.dat
		    rm -f ${OUT_FILE}

		    while read mx my mz
		    do
			EX_FILE=${EX_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
			
			${EXE_DIR}/extract_threep_h5 ${FILE} ${EX_FILE} ${tp} ${proj} ${part} ${mx} ${my} ${mz} ${CONF} ${src} ${tsink}

			cat ${EX_FILE} >> ${OUT_FILE}
		    done < ${MOM_LIST}
		done
	    done
	done
    done
done < ${SRC_LIST}
