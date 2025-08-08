#!/bin/bash

# Generar tabla CSV 'samplesinfo' a partir de directorio con fastq illumina

# Uso: ./generate_sample_table.sh /ruta/al/directorio_fastq[/]

indir="$1"

if [[ -z "$indir" ]]; then
  echo "Uso: $0 <directorio_con_fastq>" >&2
  exit 1
fi

# Encabezado
echo -e "sample,fastq1,fastq2" > samplesinfo.csv

# Buscar R1 y emparejar con R2
for r1 in "$indir"/*_R1_*.fastq.gz; do
  # Obtener nombre base sin ruta
  base=$(basename "$r1")
  
  # Extraer el ID de muestra (antes del primer "_")
  sample_id="${base%%_*}"

  # Generar nombre del R2 esperado
  r2="${r1/_R1_/_R2_}"

  # Verificar que el archivo R2 existe y completar tabla
  if [[ -f "$r2" ]]; then
    echo -e "${sample_id},${r1},${r2}" | sed 's|\/\/|\/|g' >> samplesinfo.csv
  else
    echo "⚠️  Advertencia: no se encontró R2 para $r1" >&2
  fi
done

