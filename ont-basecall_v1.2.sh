#!/bin/bash
############################################################################
# Automated basecalling for gridion runs                                   #
#                                                                          #
# Syntax: ont-basecall_v1.2.sh /path/to/run/output/folder/ barcode.csv     #
#    eg.: ont-basecall_v1.2.sh /data/220101_run1/ run1.csv                 #
#                                                                          #
# Date: 01.06.2022    Created by: Sathom                                   #
# Mod:  18.06.2022                                                         #
#                                                                          #
############################################################################

##############
### VARIABLES
##############
KIT_NAME='SQK-RBK110-96'
#FLOWCELL='FLO-MIN106'

########################
### OPTIONAL VARIABLES
########################
### default: /opt/ont/guppy/data
CFG_FILE_PATH="/opt/ont/guppy/data"
### default: dna_r9.4.1_450bps_fast.cfg
CFG_FILE_NAME="dna_r9.4.1_450bps_fast.cfg"

### type: nvidia-smi ; default: 'auto'
CUDA_DEV="cuda:all:100%"

### type: lscpu for number of threads, CPUs and sockets. THREADS= THREADS x CPU x SOCKETS
THR=12


### Switches: 1 = only basecalling; 2 = basecalling & barcoding; 3 = barcoding, basecalling and merging; 4 = all previous + renaming
SWITCH=4
A=1
B=1
C=1
D=1

##############
### PROGRAM
##############

WORKDIR=$1
BARCODE_FILE=$2

### checking working directory
STARTDIR=$(pwd)
if [   -z ${WORKDIR} ]; then echo -e "\nSyntax: $0 /path/to/output_folder/\neg. basecall.sh /data/run1_220101/"; exit 1; fi;
if [ ! -e ${WORKDIR} ]; then echo -e "\nspecified path does not exist\n"; exit 2; fi;
cd ${WORKDIR}

OPTS=''
### checking parameters
if [ -e "/opt/ont/guppy/data" ] && [[ -z ${CFG_FILE_PATH} || ! -e ${CFG_FILE_PATH} ]]; then
  CFG_FILE_PATH='/opt/ont/guppy/data';
  echo -e "\nConfiguration file path does not exists. Changed to default:\n${CFG_FILE_PATH}\n"; fi;

if [   -z ${CFG_FILE_NAME} ] && [[   -z ${KIT_NAME} &&   -z ${FLOWCELL} ]]; then echo -e "\nYou need to specify configuration file or flowcell and kit name\n"; exit 2; fi;
if [   -z ${CFG_FILE_NAME} ] && [[ ! -z ${KIT_NAME} && ! -z ${FLOWCELL} ]]; then OPTS=" --flowcell ${FLOWCELL} --kit ${KIT_NAME} "; fi;
if [ ! -z ${CFG_FILE_NAME} ] && [[ ! -z ${KIT_NAME} && ! -z ${FLOWCELL} ]]; then OPTS=" --flowcell ${FLOWCELL} --kit ${KIT_NAME} "; fi;
if [ ! -z ${CFG_FILE_NAME} ] && [[   -z ${KIT_NAME} ||   -z ${FLOWCELL} ]]; then
  if [ ! -e ${CFG_FILE_PATH} ]; then CFG_FILE_NAME='dna_r9.4.1_450bps_fast.cfg'; fi;
  OPTS=" -c ${CFG_FILE_PATH}/${CFG_FILE_NAME} ";
fi

if [ -z ${CUDA_DEV} ]; then CUDA_DEV='auto'; fi;
if [ -z ${THR} ];      then THR=2; fi;


### making output folder for FASTQ files
FAST5_PATH=`find ${WORKDIR} -type d -name 'fast5'`
OUTPUT_PATH=`echo ${FAST5_PATH} | sed -r 's/\/fast5//'`
LOG=${OUTPUT_PATH}'/_'`head -n2 ${OUTPUT_PATH}/barcode_alignment_*.tsv | tail -n1  | awk '{print $6}'`'.log' ### not always present (in not ended runs it's not)
echo "START @: "`date +%d.%m.%Y_%H:%M:%S` >>${LOG}

### checking barcode file

BARC_FILE_CHECK=(`find /data/* -type f -name ${BARCODE_FILE} 2>/dev/null `)

if [[ ${BARC_FILE_CHECK} == "" || -z ${BARC_FILE_CHECK} ]]; then echo -e "\nUnable to do barcoding, will stop before that"; ((${SWITCH}-1)); fi;

i=1;
while [[ ${i} -le ${#BARC_FILE_CHECK[@]} ]]; do
 I=$((${i}-1))
 if [[ ! -z ${BARC_FILE_CHECK[$I]} && ${BARC_FILE_CHECK[$I]} == ${OUTPUT_PATH}/${BARCODE_FILE} ]]; then chmod +x ${OUTPUT_PATH}/${BARCODE_FILE}; break; fi;
 if [[ ! -z ${BARC_FILE_CHECK[$I]} && ${BARC_FILE_CHECK[$I]} != ${OUTPUT_PATH}/${BARCODE_FILE} ]]; then cp -v ${BARC_FILE_CHECK[$I]} ${OUTPUT_PATH}/${BARCODE_FILE} && chmod +x ${OUTPUT_PATH}/${BARCODE_FILE}; break; fi;
 i=$((${i}+1))
done

FASTQ_PATH=${OUTPUT_PATH}'/Step1_fastq'

if [ ${A} -eq 1 ]; then
### STEP 1.  ### basecalling
 echo -e "Doing Step 1: Basecalling\n"
 if [ -e ${FASTQ_PATH} ]; then rm -rf ${FASTQ_PATH}; fi;
 mkdir -p -m a=rwx ${FASTQ_PATH}

 guppy_basecaller -i ${FAST5_PATH}/ -s ${FASTQ_PATH}/ --cpu_threads_per_caller ${THR} --num_callers 1 ${OPTS} --device ${CUDA_DEV} --gpu_runners_per_device ${THR} 2>&1 1>>${LOG}

 if [ ${SWITCH} -eq 1 ]; then echo -e "Only basecalling was done\n"; exit; fi;

fi

BCODES=${OUTPUT_PATH}'/Step2_barcodes'

if [ ${B} -eq 1 ]; then
 ### STEP 2.  ### barcoding
 echo -e "Doing  Step 2: Barcoding\n"
 if [ -e ${BCODES} ]; then rm -rf ${BCODES}; fi;
 mkdir -p -m a=rwx ${BCODES}

 guppy_barcoder  --input_path ${FASTQ_PATH}/pass/ --save_path ${BCODES} --config /opt/ont/guppy/data/barcoding/configuration.cfg --barcode_kits ${KIT_NAME} -t ${THR} 2>&1 1>>${LOG}

 if [ ${SWITCH} -eq 2 ]; then echo -e "Only basecalling and barcoding was done\n"; exit; fi;

fi

MERGED=${OUTPUT_PATH}'/Step3_merged'

if [ ${C} -eq 1 ]; then
 ### STEP 3.  ### merging barcodes
 echo -e "Doing Step 3: Merging barcodes\n"
 if [ -e ${MERGED} ]; then rm -rf ${MERGED}; fi;
 mkdir -p -m a=rwx ${MERGED}
 
 for i in $(ls -l ${BCODES} | grep ^d | awk '$9 ~/^barcode/ {print $9}');
   do cat ${BCODES}/${i}/*.fastq >> ${MERGED}/${i}.fastq;
 done

 if [ ${SWITCH} -eq 3 ]; then echo -e "Only basecalling, barcoding and barcodes merge was done\n"; exit; fi;
fi

RENAMED=${OUTPUT_PATH}'/Step4_renamed'

if [ ${D} -eq 1 ]; then
 ### STEP 4. ### decoding barcodes
 echo -e "Doing Step 4: Decoding barcodes\n"
 if [ -e ${RENAMED} ]; then rm -rf ${RENAMED}; fi;
 mkdir -p -m a=rwx ${RENAMED}

 ### reading barcode file & renaming files
 for L in $(cat ${OUTPUT_PATH}/${BARCODE_FILE}); do
  Z=''
  K=`echo $L | cut -d ',' -f1`
  V=`echo $L | cut -d ',' -f2`
  if [ ${K} -lt 10 ]; then Z='0'; fi;
  cp ${MERGED}/barcode${Z}${K}.fastq ${RENAMED}/${V}_ONT.fastq
 done

 cd ${RENAMED}
 parallel --version 2>/dev/null
 if [ $? -eq 0 ]; then ls *.fastq | parallel -j ${THR} -I% "gzip -9 %"; fi;
 gzip -9 *.fastq

 if [ ${SWITCH} -eq 4 ]; then echo -e "All steps done\n"; fi;
fi

cd ${STARTDIR}
##########################################################
### has to be added to all premature ends!!!
echo "END @: "`date +%d.%m.%Y_%H:%M:%S` >>${LOG}

exit 0
