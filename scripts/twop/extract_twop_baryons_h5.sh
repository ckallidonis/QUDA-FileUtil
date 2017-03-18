#!/bin/bash

if [ $# -ne 10 ]
then
    echo "Usage: $0 <h5_file-prefix> <out_dir> <src_list> <mom_list> <conf> <T> <Qsq> <exe_dir> <HighMomForm:1 LowMomForm:0> <baryon, see list below>"
    printf "  0: all\n  1: nucl_nucl\n  2: nucl_roper\n  3: roper_nucl\n  4: roper_roper\n  5: deltapp_deltamm_11\n  6: deltapp_deltamm_22\n  7: deltapp_deltamm_33\n  8: deltap_deltaz_11\n  9: deltap_deltaz_22\n 10: deltap_deltaz_33\n"
    exit
fi

H5_PRE=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
T=$6
Qsq=$7
EXE_DIR=$8
HighMomForm=${9}
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


if [ $HighMomForm -eq 0 ]
then

    while read x y z t
    do 
	FILE=${H5_PRE}${CONF}_baryons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}
	EX_DIR=${OUT_DIR}/${src}
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

	rm -rf ${EX_DIR}

    done < ${SRC_LIST}

fi

if [ $HighMomForm -eq 1 ]
then

    while read x y z t
    do 
	FILE=${H5_PRE}${CONF}_baryons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}
	EX_DIR=${OUT_DIR}/${src}
	mkdir -p ${EX_DIR}

	OUT_FILE=${OUT_DIR}/twop.${CONF}.baryons.SS.${x}.${y}.${z}.${t}.dat

	rm -f ${OUT_FILE}

	for tp in ${BARYON_LIST}
	do
	    EX_FILE=${EX_DIR}/twop.${CONF}.baryons.${tp}.SS.${x}.${y}.${z}.${t}.dat
		
	    ${EXE_DIR}/extract_twop_baryons_HighMomForm_h5 ${FILE} ${EX_FILE} ${tp} ${CONF} ${src} ${T} 0

	    cat ${EX_FILE} >> ${OUT_FILE}
	done

	rm -rf ${EX_DIR}

    done < ${SRC_LIST}

fi
