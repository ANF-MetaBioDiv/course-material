#!/bin/bash/

FORWARD_IN=$1
REVERSE_IN=$2
PRIMER_F=$3
PRIMER_R=$4
OUT_DIR=$5
MIN_LENGTH=$6
LOG_DIR=$7

FORWARD_OUT="${OUT_DIR}/$(basename ${FORWARD_IN})"
REVERSE_OUT="${OUT_DIR}/$(basename ${REVERSE_IN})"

SAMPLE=$(basename ${FORWARD_IN/_R[12]*.fastq.gz/})

LOG="${LOG_DIR}/${SAMPLE}_primer_trimming.log"
TMP_LOG=$(mktemp --tmpdir=".")

cutadapt \
	-g ${PRIMER_F} \
	-G ${PRIMER_R} \
	--report=minimal \
	--discard-untrimmed \
	--minimum-length ${MIN_LENGTH} \
	--no-indels \
	-o ${FORWARD_OUT} \
	-p ${REVERSE_OUT} \
	${FORWARD_IN} \
	${REVERSE_IN} 1> ${TMP_LOG}

awk -v a="${SAMPLE}" \
	'BEGIN {OFS="\t"}; NR==2{print a,$0}' \
	${TMP_LOG} > ${LOG}

rm -f ${TMP_LOG}
