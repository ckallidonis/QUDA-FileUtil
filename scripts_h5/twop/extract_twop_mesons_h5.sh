#!/bin/bash

if [ $# -ne 7 ]
then
    echo "Usage: $0 <extract_dir> <src_list> <mom_list> <conf> <T> <Qsq> <exe_dir>"
    exit
fi

DIR=$1
SRC_LIST=$2
MOM_LIST=$3
CONF=$4
T=$5
Qsq=$6
EXE_DIR=$7

while read x y z t
do 
    FILE=${DIR}/twop.${CONF}_mesons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}

    OUT_DIR=${DIR}/${src}

    mkdir -p ${OUT_DIR}

    for tp in pseudoscalar scalar g5g1 g5g2 g5g3 g5g4 g1 g2 g3 g4
    do
	while read mx my mz
	do
	    OUT_FILE=${OUT_DIR}/twop.${CONF}.mesons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
	    
	    ${EXE_DIR}/extract_twop_mesons_h5 ${FILE} ${OUT_FILE} ${tp} ${mx} ${my} ${mz} ${CONF} ${src} ${T} 0
	done < ${MOM_LIST}
    done

done < ${SRC_LIST}
