#! /bin/bash

## 1) Obtener dataset asignado por el pipeline (para las muestras en que haya):
RUN=$1
inpath="$PWD/results/${RUN}/nextclade_variables/"
outpath="$PWD/work/"
ls ${inpath} | while read folder ; do echo ${folder} | tr "\n" "\t" ; ls ${inpath}${folder} ; done > ${outpath}sample_dataset.tmp.tsv

## 2) run nextclade
# Nextclade
path_fasta="$PWD/results/${RUN}/hemagglutinin/"
path_dataset="$PWD/results/${RUN}/nextclade_datasetget/"
adate=$(date '+%y%m%d') 

cat ${outpath}sample_dataset.tmp.tsv \
	| while read sample dataset; do 
		nextclade run --input-dataset ${path_dataset}${dataset}* \
					  --output-all ${output}${adate}_Nextclade_${RUN} \
					  --output-basename ${sample} ${path_fasta}${sample}_HA.fasta
done


