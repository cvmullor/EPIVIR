# EPIVIR

*Pipeline* para la vigilancia genómica de las infecciones respiratorias agudas (IRAs) en la Comunitat Valenciana, diseñado para el análisis de diferentes organismos (gripe, VRS) y tecnologías de secuenciación (Illumina, Oxford Nanopore)

**EPIVIR** es una adaptación del *pipeline* **walkercreek**. Para más información sobre este y las modificaciones que se han implementado, consultar la sección [About this fork](https://github.com/cvmullor/EPIVIR/tree/master?tab=readme-ov-file#-about-this-fork).

---

## Puesta a punto EPIVIR
### Instalación de Conda
```bash
[pending]
```

### Descargar repositorio
```bash
git clone https://github.com/cvmullor/EPIVIR.git
```

### Instalación de dependencias
Se instalarán las depencias necesarias para utilizar EPIVIR en un entorno de Conda llamado `epivir-env`.
```bash
conda env create -y -f EPIVIR/env.yaml
```

## Lanzar EPIVIR
Activar el entorno de Conda.
```bash
conda activate epivir-env
```

Generar tabla *samplesinfo.csv* que EPIVIR emplea como input, indicando la ruta al directorio con los archivos fastq de las muestras y la plataforma (illumina|nanopore) de secuenciación. La tabla se generará en el directorio de trabajo actual.
```bash
./EPIVIR/create_samplesinfo.sh <directorio_fastq> <plataforma>
```

Ejecutar script principal, indicando organismo (FLU|RSV), plataforma de secuenciación (illumina|nanopore), nombre de la carrera y ruta a la tabla *samplesinfo.csv*.
```bash
./EPIVIR/run_epivir.sh <organismo> <plataforma> <carrera> <samplesinfo>
```

Ejemplo:
```bash
./EPIVIR/run_epivir.sh FLU illumina 250101_Test_FLU samplesinfo.csv
```

### Resultados
+ Directorio `results/<carrera>/`: se encuentran los ficheros de resultados de cada análisis independiente realizado por EPIVIR.
+ Tabla `results/<carrera>/<fecha>_RESULTS_<carrera>.csv`: resumen del análisis por muestra. Se incluyen campos correspondientes a los metadatos del aislado o paciente, a rellenar por el hospital.

---


## 🔁 About this fork

This repository is a fork of [**walkercreek**](https://github.com/UPHL-BioNGS/walkercreek), developed by Tom Iverson ([@tives82](https://github.com/tives82)) and contributors.

*walkercreek* is a Nextflow pipeline designed to accommodate different sequencing technologies (Illumina or Nanopore), sample types (clinical or wastewater), and viral targets (Influenza or RSV).

This fork has been adapted for **genomic surveillance of acute respiratory infections (ARIs) in the Comunitat Valenciana, Spain**.

The original project is licensed under the [MIT License](https://opensource.org/licenses/MIT).  
Please see the [`LICENSE`](LICENSE) file for full licensing details.
