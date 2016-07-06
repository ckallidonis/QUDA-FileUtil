#!/bin/bash

if [ $# -ne 8 ]
then
    echo "Usage: $0 <extract_dir> <out_dir> <src_list> <mom_list> <conf> <T> <Qsq> <exe_dir>"
    exit
fi

DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
T=$6
Qsq=$7
EXE_DIR=$8

while read x y z t
do 
    FILE=${DIR}/twop.${CONF}_mesons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}
    EX_DIR=${DIR}/${src}
    mkdir -p ${EX_DIR}

    OUT_FILE=${OUT_DIR}/twop.${CONF}.mesons.SS.${x}.${y}.${z}.${t}.dat

    rm -f ${OUT_FILE}

    for tp in pseudoscalar scalar g5g1 g5g2 g5g3 g5g4 g1 g2 g3 g4
    do
	while read mx my mz
	do
	    EX_FILE=${EX_DIR}/twop.${CONF}.mesons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
	    
	    ${EXE_DIR}/extract_twop_mesons_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${CONF} ${src} ${T} 0

            cat ${EX_FILE} >> ${OUT_FILE}

	done < ${MOM_LIST}
    done

done < ${SRC_LIST}
