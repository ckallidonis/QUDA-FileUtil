#!/bin/bash

if [ $# -ne 10 ]
then
    echo "Usage: $0 <h5_file-prefix> <out_dir> <src_list> <mom_list> <conf> <T> <Qsq> <exe_dir> <HighMomForm:1 LowMomForm:0> <meson, see list below"
    printf "  0: all\n  1: pseudoscalar\n  2: scalar\n  3: g5g1\n  4: g5g2\n  5: g5g3\n  6: g5g4\n  7: g1\n  8: g2\n  9: g3\n 10: g4\n"
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
mes=${10}

LIST[1]="pseudoscalar"
LIST[2]="scalar"
LIST[3]="g5g1"
LIST[4]="g5g2"
LIST[5]="g5g3"
LIST[6]="g5g4"
LIST[7]="g1"
LIST[8]="g2"
LIST[9]="g3"
LIST[10]="g4"

if [ $mes -eq 0 ]
then
    MESON_LIST="pseudoscalar scalar g5g1 g5g2 g5g3 g5g4 g1 g2 g3 g4"
fi

if [ $mes -gt 0 ]
then
    MESON_LIST=${LIST[$mes]}
fi


if [ ${HighMomForm} -eq 0 ]
then

    while read x y z t
    do 
	FILE=${H5_PRE}${CONF}_mesons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}
	EX_DIR=${OUT_DIR}/${src}
	mkdir -p ${EX_DIR}

	OUT_FILE=${OUT_DIR}/twop.${CONF}.mesons.SS.${x}.${y}.${z}.${t}.dat

	rm -f ${OUT_FILE}

	for tp in ${MESON_LIST}
	do
	    while read mx my mz
	    do
		EX_FILE=${EX_DIR}/twop.${CONF}.mesons.${tp}.${mx}_${my}_${mz}.SS.${x}.${y}.${z}.${t}.dat
		
		${EXE_DIR}/extract_twop_mesons_h5 ${FILE} ${EX_FILE} ${tp} ${mx} ${my} ${mz} ${CONF} ${src} ${T} 0

		cat ${EX_FILE} >> ${OUT_FILE}

	    done < ${MOM_LIST}
	done

	rm -rf ${EX_DIR}

    done < ${SRC_LIST}

fi


if [ ${HighMomForm} -eq 1 ]
then

    while read x y z t
    do 
	FILE=${H5_PRE}${CONF}_mesons_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}
	EX_DIR=${OUT_DIR}/${src}
	mkdir -p ${EX_DIR}

	OUT_FILE=${OUT_DIR}/twop.${CONF}.mesons.SS.${x}.${y}.${z}.${t}.dat

	rm -f ${OUT_FILE}

	for tp in ${MESON_LIST}
	do
	    EX_FILE=${EX_DIR}/twop.${CONF}.mesons.${tp}.SS.${x}.${y}.${z}.${t}.dat
	    
	    ${EXE_DIR}/extract_twop_mesons_HighMomForm_h5 ${FILE} ${EX_FILE} ${tp} ${CONF} ${src} ${T} 0
	    
	    cat ${EX_FILE} >> ${OUT_FILE}
	done

	rm -rf ${EX_DIR}

    done < ${SRC_LIST}

fi
