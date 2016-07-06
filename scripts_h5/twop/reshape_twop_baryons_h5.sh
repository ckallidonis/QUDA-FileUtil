#!/bin/bash

if [ $# -ne 5 ]
then
    echo "Usage: $0 <extract_dir> <out_dir> <src_list> <mom_list> <conf>"
    exit
fi

DIR=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5

while read x y z t
do
    src=sx${x}sy${y}sz${z}st${t}
    EX_DIR=${DIR}/${src}

    OUT_FILE=${OUT_DIR}/twop.${CONF}.baryons.SS.${x}.${y}.${z}.${t}.dat	

    rm -f ${OUT_FILE}

    for tp in nucl_nucl nucl_roper roper_nucl roper_roper deltapp_deltamm_11 deltapp_deltamm_22 deltapp_deltamm_33 deltap_deltaz_11 deltap_deltaz_22 deltap_deltaz_33
    do
	
	while read mx my mz
	do
	    FILE=${EX_DIR}/twop.${CONF}.baryons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
	          	    
	    cat ${FILE} >> ${OUT_FILE}
	    
	done < ${MOM_LIST}
    done

    echo $src

done < ${SRC_LIST}
