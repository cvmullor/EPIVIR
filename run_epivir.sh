#! /bin/bash

# Argumentos
epidir="$PWD/EPIVIR"
basedir="$PWD/results"
workdir="$PWD/work"

virus=$1    # FLU o RSV
platform=$2 # illumina o nanopore
run=$3		# nombre del run
ssinfo=$4   # ruta a la SamplesInfo.csv

run_epivir=true

if [[ "$virus" == "FLU" && "$platform" == "illumina" ]]; then
	mode="flu_illumina"

elif [[ "$virus" == "RSV" && "$platform" == "illumina" ]]; then
	mode="rsv_illumina"

elif [[ "$virus" == "FLU" && "$platform" == "nanopore" ]]; then
	mode="flu_nanopore"

elif [[ "$virus" == "RSV" && "$platform" == "nanopore" ]]; then
	mode="rsv_nanopore"

else
	run_epivir=false
	echo "ERROR Argumento 'virus' o 'platform' no válidos"
fi


### Análisis EPIVIR ###

# Previamente, cargar entorno: "conda activate epivir-env"

# Ejecutar pipeline:
if [[ "$run_epivir" == "true" ]]; then 

	nextflow run EPIVIR -cache true -profile singularity \
		--platform ${mode} \
		--input ${ssinfo} \
		--outdir ${basedir}/${run}

	# Generar reporte
	if [[ "$platform" == "illumina" ]]; then

		Rscript ${epidir}/custom_report.R ${virus} ${run} ${basedir} ${epidir} ${workdir}

	fi
fi

