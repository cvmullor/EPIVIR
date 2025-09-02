#!/bin/bash

# Generar tabla CSV 'samplesinfo' a partir de directorio con fastq
# Compatible con Illumina (PE) y Nanopore
#
# Uso:
#   ./generate_sample_table.sh <directorio_fastq> <plataforma>
#   plataforma = illumina | nanopore

indir="$1"
platform="$2"

if [[ -z "$indir" || -z "$platform" ]]; then
  echo "Uso: $0 <directorio_con_fastq> <illumina|nanopore>" >&2
  exit 1
fi

# Normalizar ruta (quitar // dobles)
indir=$(echo "$indir" | sed 's|//*|/|g')

# Illumina
if [[ "$platform" == "illumina" ]]; then
  echo -e "sample,fastq1,fastq2" > samplesinfo.csv

  for r1 in "$indir"/*_R1_*.fastq.gz; do
    # Evitar glob vacío
    [[ -e "$r1" ]] || continue
    
    base=$(basename "$r1")
    sample_id="${base%%_*}"
    r2="${r1/_R1_/_R2_}"

    if [[ -f "$r2" ]]; then
      echo -e "${sample_id},${r1},${r2}" | sed 's|//*|/|g' >> samplesinfo.csv
    else
      echo "⚠️  Advertencia: no se encontró R2 para $r1" >&2
    fi
  done

# Nanopore
elif [[ "$platform" == "nanopore" ]]; then
  echo -e "sample,fastq" > samplesinfo.csv

  for fq in "$indir"/*.fastq.gz; do
    [[ -e "$fq" ]] || continue
    
    base=$(basename "$fq")
    sample_id="${base%%.*}"
    echo -e "${sample_id},${fq}" | sed 's|//*|/|g' >> samplesinfo.csv
  done

else
  echo "Plataforma no reconocida: $platform. Usa 'illumina' o 'nanopore'." >&2
  exit 1
fi

echo "Tabla generada: samplesinfo.csv"

