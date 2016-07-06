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
    FILE=${DIR}/twop.${CONF}_baryons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}

    OUT_DIR=${DIR}/${src}

    mkdir -p ${OUT_DIR}

    for tp in nucl_nucl nucl_roper roper_nucl roper_roper deltapp_deltamm_11 deltapp_deltamm_22 deltapp_deltamm_33 deltap_deltaz_11 deltap_deltaz_22 deltap_deltaz_33
    do
	while read mx my mz
	do
	    OUT_FILE=${OUT_DIR}/twop.${CONF}.baryons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
	    
	    ${EXE_DIR}/extract_twop_baryons_h5 ${FILE} ${OUT_FILE} ${tp} ${mx} ${my} ${mz} ${CONF} ${src} ${T} 0
	done < ${MOM_LIST}
    done

done < ${SRC_LIST}
