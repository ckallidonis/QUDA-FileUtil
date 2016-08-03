#!/bin/bash

if [ $# -ne 10 ]
then
    echo "Usage: $0 <extract_dir> <h5_file-prefix> <out_dir> <src_list> <mom_list> <conf> <T> <Qsq> <exe_dir>  <baryon, see list below>"
    printf "  0: all\n  1: nucl_nucl\n  2: nucl_roper\n  3: roper_nucl\n  4: roper_roper\n  5: deltapp_deltamm_11\n  6: deltapp_deltamm_22\n  7: deltapp_deltamm_33\n  8: deltap_deltaz_11\n  9: deltap_deltaz_22\n 10: deltap_deltaz_33\n"
    exit
fi

DIR=$1
H5_PRE=$2
OUT_DIR=$3
SRC_LIST=$4
MOM_LIST=$5
CONF=$6
T=$7
Qsq=$8
EXE_DIR=$9
bar=${10}

LIST[1]="nucl_nucl"
LIST[2]="nucl_roper"
LIST[3]="roper_nucl"
LIST[4]="roper_roper"
LIST[5]="deltapp_deltamm_11"
LIST[6]="deltapp_deltamm_22"
LIST[7]="deltapp_deltamm_33"
LIST[8]="deltap_deltaz_11"
LIST[9]="deltap_deltaz_22"
LIST[10]="deltap_deltaz_33"

if [ $bar -eq 0 ]
then
    BARYON_LIST="nucl_nucl nucl_roper roper_nucl roper_roper deltapp_deltamm_11 deltapp_deltamm_22 deltapp_deltamm_33 deltap_deltaz_11 deltap_deltaz_22 deltap_deltaz_33"
fi

if [ $bar -gt 0 ]
then
    BARYON_LIST=${LIST[$bar]}
fi

while read x y z t
do 
    FILE=${DIR}/${H5_PRE}.${CONF}_baryons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

    src=sx${x}sy${y}sz${z}st${t}
    EX_DIR=${DIR}/${src}
    mkdir -p ${EX_DIR}

    OUT_FILE=${OUT_DIR}/twop.${CONF}.baryons.SS.${x}.${y}.${z}.${t}.dat

    rm -f ${OUT_FILE}

    for tp in ${BARYON_LIST}
    do
	while read mx my mz
	do
	    EX_FILE=${EX_DIR}/twop.${CONF}.baryons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
	    
	    ${EXE_DIR}/extract_twop_baryons_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${CONF} ${src} ${T} 0

	    cat ${EX_FILE} >> ${OUT_FILE}

	done < ${MOM_LIST}
    done

done < ${SRC_LIST}
