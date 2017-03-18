#!/bin/bash

if [ $# -ne 11 ]
then
    echo "Usage: $0 <h5_prefix> <out_dir> <src_list> <mom_list> <conf> <tsink> <proj> <Qsq> <exe_dir> <HighMomForm:1 LowMomForm:0> <type, see below>"
    printf " 0: all\n 1: ultra_local\n 2: noether\n 3: oneD\n"
    exit
fi

PREFIX=$1
OUT_DIR=$2
SRC_LIST=$3
MOM_LIST=$4
CONF=$5
TSINK=$6
PROJ=$7
Qsq=$8
EXE_DIR=$9
HighMomForm=${10}
EX_TYPE=${11}

if [ ${EX_TYPE} -eq 0 ]
then
    EX_LIST="ultra_local noether oneD"
fi
if [ ${EX_TYPE} -eq 1 ]
then
    EX_LIST="ultra_local"
fi
if [ ${EX_TYPE} -eq 2 ]
then
    EX_LIST="noether"
fi
if [ ${EX_TYPE} -eq 3 ]
then
    EX_LIST="oneD"
fi


if [ ${HighMomForm} -eq 0 ]
then

    while read x y z t
    do 
	FILE=${PREFIX}${CONF}_neutron_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}

	EX_DIR=${OUT_DIR}/${src}

	mkdir -p ${EX_DIR}

	for tsink in $TSINK
	do 
	    for proj in $PROJ
	    do
		for part in up down
		do
		    for tp in ${EX_LIST}
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

	rm -rf ${EX_DIR}

    done < ${SRC_LIST}

fi



if [ ${HighMomForm} -eq 1 ]
then

    while read x y z t
    do 
	FILE=${PREFIX}${CONF}_neutron_Qsq${Qsq}_SS.${x}.${y}.${z}.${t}.h5

	src=sx${x}sy${y}sz${z}st${t}

	for tsink in $TSINK
	do 
	    for proj in $PROJ
	    do
		for part in up down
		do
		    for tp in ${EX_LIST}
		    do
			OUT_FILE=${OUT_DIR}/threep.${CONF}_tsink${tsink}_proj${proj}.neutron.${part}.${tp}.SS.${x}.${y}.${z}.${t}.dat

			${EXE_DIR}/extract_threep_HighMomForm_h5 ${FILE} ${OUT_FILE} ${tp} ${proj} ${part} ${CONF} ${src} ${tsink}
		    done
		done
	    done
	done

    done < ${SRC_LIST}

fi
