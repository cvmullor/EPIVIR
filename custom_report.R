## Generar report resultados walkercreek para FLU/RSV

## ARGUMENTS

# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments must be supplied (virus, runid)", call.=FALSE)
} else {
  virus    <- args[1]  # FLU / RSV
  RUN      <- args[2]  # e.g.: "250203_EPIM197_RSV"
  basedir0 <- args[3]  # directorio con los resultados de EPIVIR
  epidir   <- args[4]  # directorio EPIVIR
  workdir  <- args[5]  # directorio "work"
}

## PATH
basedir <- paste0(basedir0, "/", RUN, "/")
outdir0 <- paste0(basedir0, "/", RUN, "/")

## GENERAR REPORT ##

############# RSV #################

if ( virus == "RSV") {
    
  ###  1) read summary report, select columns ###
  
  indir <- "SUMMARY_REPORT/"
  afile <- "summary_report.tsv"
  
  f1 <- read.table(paste0(basedir, indir, afile), sep="\t", header = T) # read table summary_report.tsv
  
  # Seleccionar columnas
  df1 <- f1[, c('Sample', 'IRMA_consensus_N_count', 'reads_before_trimming', 'reads_after_trimming',
                'kraken2.Respiratory.Syncityal.Virus.percentage', 'kraken2.Homo.sapiens.percentage',
                'IRMA_subtype', 'clade', 'Nextclade_qc.overallScore', 'Nextclade_qc.overallStatus')]
  
  df1[df1==""] <- NA # replace blank space with NA
  
  df1 <- df1[order(df1$Sample),] # Ordenar por Sample (Hospital.sample.number)
  
  
  ### 2) coverage, depth y mapped reads ###
  
  indir <- "reports/"
  afile <- "merged_bam_coverage_results.tsv"
  
  f2 <- read.table(paste0(basedir, indir, afile), sep="\t", header = T, row.names = 1) # read table
  
  # !! --> Puede haber menos muestras analizadas que el report final.
  
  # (como solo hay un segment, a diferencia de FLU): seleccionar directamente las columnas.
  df2 <- cbind(rownames(f2), f2[, c(2, 3, dim(f2)[2])])
  
  colnames(df2) <- c("Sample","Mapped_Reads","Mean_Depth","Coverage")
  
  df2 <- df2[order(df2$Sample),] # Ordenar por Sample (Hospital.sample.number)
  

  
  ### 3) Extraer resultados Nextclade que no se incluyen en Summary Report ###
 
  ## Leer resultados nextclade por muestra:
  
  indir2 <- "nextclade_run"
  nc_list <- list.files(paste0(basedir, indir2), pattern=".tsv", recursive = T, full.names = T)
  
  # Cada muestra pueden tener num. cols diferente: extraemos en bucle las columnas de cada tabla individual.
  
  df3 <- t(as.data.frame(rep(NA, 4)))
  colnames(df3) <- c('seqName', 'totalFrameShifts', 'frameShifts', 'aaInsertions')
  rownames(df3) <- c('0')
  
  for (ti in nc_list) {
    fi <- read.table(ti, sep="\t", header = T) # read table
    dfi <- fi[, c('seqName', 'totalFrameShifts', 'frameShifts', 'aaInsertions')] # select cols
    
    # (Para RSV) Sample name no está en la tabla: extraemos del nombre del archivo
    sepi <- unlist(strsplit(ti, "/"))
    filei <- sepi[length(sepi)]
    namei <- unlist(strsplit(filei, ".", fixed = T))[1]
    dfi$seqName <- namei
    
    df3 <- rbind(df3, dfi)
  }
  
  df3 <- df3[-1,] # eliminar primera linea en blanco
  
  colnames(df3) <- c('Sample', 'totalFrameShifts', 'frameShifts', 'aaInsertions')
  
  df3 <- df3[order(df3$Sample),]  # Ordenar por Sample (Hospital.sample.number)
  
  
  ### 4) Merge results ###
  dim(df1)
  dim(df2)
  dim(df3)
  
  # Merge tables by sample name
  df4.0 <- merge(df1, df2, by.x="Sample", by.y="Sample", all = TRUE)
  df4 <- merge(df4.0, df3, by.x="Sample", by.y="Sample", all = TRUE)
  
  df4[df4==""] <- NA # replace blank space with NA
  
  
  ### 5) Generar columnas vacías / con valores constantes ###
  nr <- nrow(df4)
  
  ID <- rep("", nr)
  SIP <- rep("", nr)
  Hospital <- rep("", nr)       # completar
  Reception.date <- rep("", nr) # completar
  Hospital.date <- rep("", nr)  # completar
  Relevance <- rep("VIGILANCIA", nr)
  Comments <- rep("", nr)
  Gender <- rep("", nr)
  Age <- rep("", nr)
  Status <- rep("", nr) # Completar
  runid <- rep(RUN, nr)
  
  trimmed_reads <- as.numeric(unlist(lapply(strsplit(df4$reads_after_trimming, " "), '[', 1))) # Separete number from %
  percent_missed_reads <- 100 - round(( as.numeric(df4$Mapped_Reads) / trimmed_reads)*100, 2)
  
  version <- rep("Nextclade v3.13.0", nr)
  GISAID <- rep("", nr)
  
  
  ### 6) Ordenar columnas
  
  df5 <- as.data.frame(cbind(
               ID, df4$Sample, SIP, Hospital, Reception.date, Hospital.date,
               Relevance, Comments, Gender, Age, Status, df4$Coverage,
               df4$Mean_Depth, runid, df4$IRMA_consensus_N_count, df4$reads_before_trimming,
               trimmed_reads, df4$Mapped_Reads, percent_missed_reads,
               df4$kraken2.Respiratory.Syncityal.Virus.percentage,
               df4$kraken2.Homo.sapiens.percentage, df4$IRMA_subtype, df4$clade, version, 
               df4$Nextclade_qc.overallScore, df4$Nextclade_qc.overallStatus,
               df4$totalFrameShifts, df4$frameShifts, df4$aaInsertions, GISAID))
  
  
  
  colnames(df5) <- c("ID","Hospital.sample.number","SIP","Hospital","Reception.date",
                     "Hospital.date","Relevance","Comments","Gender","Age","Status",
                     "Coverage","Mean_Depth","runid","N_count","reads",
                     "reads_after_trimming","mapped_reads","percent_missed_reads",
                     "percent_Respiratory.Syncityal.Virus","percent_Homo.sapiens",
                     "type","clade","version","qc.overallScore","qc.overallStatus",
                     "totalFrameShifts","frameShifts","aaInsertions","GISAID")
  
  
  ## Calcular Status de cada muestra
  
  df5$N_count[is.na(df5$N_count)] <- 15277   # Si N count = NA, se introduce el tamaño de la referencia
  df5$Mean_Depth[is.na(df5$Mean_Depth)] <- 0
  df5$Coverage[is.na(df5$Coverage)] <- 0
  df5$Status[as.numeric(df5$Coverage) < 90 | df5$N_count > 3000] <- "Failed_QC"
  df5$Status[as.numeric(df5$Coverage) > 90 & df5$N_count < 3000] <- "Sequenced"
  
  ### 8) Write table
  di <- format(Sys.time(), "%y%m%d") # date of analysis
  
  outfile = paste0(outdir0, di, "_RESULTS_", RUN, ".csv")
  write.table(df5, outfile, row.names = F, sep=";")


  
  
############# FLU #################
  
} else if (virus == "FLU") {
  
  ###  1) read summary report, select columns ###
  
  indir <- "SUMMARY_REPORT/"
  afile <- "summary_report.tsv"
  
  f1 <- read.table(paste0(basedir, indir, afile), sep="\t", header = T) # read table summary_report.tsv
  
  # Extraer columnas de summary report
  df1 <- f1[, c('Sample', 'IRMA_consensus_N_count', 'reads_before_trimming', 'reads_after_trimming',
                'kraken2.Influenza.A.percentage', 'kraken2.Influenza.B.percentage', 'kraken2.Homo.sapiens.percentage',
                'IRMA_type', 'abricate_InsaFlu_type', 'IRMA_subtype', 'abricate_InsaFlu_subtype',
                'clade', 'subclade', 'Nextclade_qc.overallScore', 'Nextclade_qc.overallStatus')]
  
  df1[df1==""] <- NA # replace blank space with NA
  
  df1 <- df1[order(df1$Sample),] # Ordenar por Sample
  
  
  ### 2) coverage, depth y mapped reads ###
  
  indir <- "reports/"
  afile <- "merged_bam_coverage_results.tsv"
  
  f2 <- read.table(paste0(basedir, indir, afile), sep="\t", header = T, row.names = 1) # read table
  
  # !! --> No hay resultado de 1 negativo!!
  
  # Separar resultados coverage de mapped+depth
  f2.mapped_dep <- f2[, 1:24]
  f2.coverage   <- f2[, 25:dim(f2)[2]]
  
  # Get segment names
  segment_names <- c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2")
  df2 <- as.data.frame(rownames(f2))
  df2[,2:9] <- NA # create 8 empty columnes
  colnames(df2) <- c("Sample", segment_names)
  rownames(df2) <- df2$Sample
  
  # Tables for depth, coverage
  df2.cov <- df2
  df2.dep <- df2
  df2.map <- df2
  
  ### Loop ids and segments
  for ( id in rownames(f2) ) {
    for (iseg in segment_names){
      
      # reads mapped, mean depth
      idx <- which(grepl(paste0("_",iseg), f2.mapped_dep[id,])) # grep idx 1 segment
      
      if (length(idx)==0){
        df2.dep[id, iseg] <- NA # Si no hay valores, introducir NA
        df2.map[id, iseg] <- NA 
        
      } else {
        value_dm <- f2.mapped_dep[id, c(idx:(idx+2))]  # obtener los valores asociados de depth y mapped reads
        
        df2.dep[id, iseg] <- value_dm[3] # introducir valor depth por muestra y segmento
        df2.map[id, iseg] <- value_dm[2] # introducir valor n_reads_mapeadas por muestra y segmento
      }
      
      
      # Coverage
      idx <- which(grepl(paste0("_",iseg), f2.coverage[id,])) # grep idx 1 segment
      
      if (length(idx)==0){
        df2.cov[id, iseg] <- NA # Si no hay valores, introducir NA
        
      }  else {
        value_co <- f2.coverage[id, c(idx,(idx+3))]  # obtener los valores asociados a coverage
        
        df2.cov[id, iseg] <- value_co[2] # introducir valor percent_coverage
        
      }
    }
  }
  
  
  # Def dataframes to merge
  colnames(df2.cov) <- c("Sample", paste0(segment_names, ".Coverage"))   # Coverage
  colnames(df2.dep) <- c("Sample", paste0(segment_names, ".Mean_Depth")) # Mean Depth
  colnames(df2.map) <- c("Sample", paste0(segment_names, ".Mapped_Reads")) # Mean Depth
  
  # Total Mapped reads
  df2.map[is.na(df2.map)] <- 0
  df2.total_map <- as.data.frame(cbind(df2.map[,1], rowSums(df2.map[,2:8])))
  colnames(df2.total_map) <- c("Sample", "mapped_reads")
  
  
  # Sort by Sample
  df2.cov <- df2.cov[order(df2.cov$Sample),]
  df2.dep <- df2.dep[order(df2.dep$Sample),]
  df2.total_map <- df2.total_map[order(df2.total_map$Sample),]
  
  
  ### 3) Extraer resultados Nextclade no incluidos en Summary Report ###
  
  
  ####>>> Relanzar NC, parsear <<<####
  
  # run nextclade
  extra_dir <- paste0(workdir, "/")
  com1 <- paste0(epidir, "/run_Nextclade.sh ", RUN)
  system(command = com1)
  
  # mover resultsods
  di <- format(Sys.time(), "%y%m%d") # date of analysis
  #outrun <- paste0(extra_dir, di, "_Nextclade_", RUN, ".csv")
  outrun <- paste0(di, "_Nextclade_", RUN)
  
  com2 <- paste0("mv ", outrun, " ", extra_dir)
  system(command = com2)
  
  
  # parsear salida nextclade extra
  indir2 <- paste0(extra_dir, outrun)
  nc_list <- list.files(paste0(extra_dir, outrun), pattern=".tsv", recursive = T, full.names = T)
  
  # Cada muestra pueden tener num. cols diferente: extraemos en bucle las columnas de cada tabla individual.
  df3 <- t(as.data.frame(rep(NA, 8)))
  colnames(df3) <- c('seqName', 'clade', 'subclade', 'qc.overallScore', 'qc.overallStatus', 'totalFrameShifts', 'frameShifts', 'aaInsertions')
  rownames(df3) <- c('0')
  
  for (ti in nc_list) {
    
    fi <- read.table(ti, sep="\t", header = T) # read table
    
    if ( "subclade" %in% colnames(fi) ) {
      
      dfi <- fi[, c('seqName', 'clade', 'subclade', 'qc.overallScore', 'qc.overallStatus', 'totalFrameShifts', 'frameShifts', 'aaInsertions')] # select cols
      
      # (Para RSV) Sample name no está en la tabla: extraemos del nombre del archivo
      sepi <- unlist(strsplit(ti, "/"))
      filei <- sepi[length(sepi)]
      namei <- unlist(strsplit(filei, ".", fixed = T))[1]
      dfi$seqName <- namei
      
      df3 <- rbind(df3, dfi)
      
    } else if ( ! "subclade" %in% colnames(fi)) {
      
      dfi_0 <- fi[, c('seqName', 'clade', 'qc.overallScore', 'qc.overallStatus', 'totalFrameShifts', 'frameShifts', 'aaInsertions')] # select cols
      subclade <- c(NA)
      dfi <- cbind(dfi_0[,1:2], subclade, dfi_0[,3:7])
      
      # (Para RSV) Sample name no está en la tabla: extraemos del nombre del archivo
      sepi <- unlist(strsplit(ti, "/"))
      filei <- sepi[length(sepi)]
      namei <- unlist(strsplit(filei, ".", fixed = T))[1]
      dfi$seqName <- namei
      
      df3 <- rbind(df3, dfi)
      
    }
  }  
    

  # df3 <- df3[-1,] # eliminar primera linea en blanco
  
  colnames(df3) <- c('Sample', 'clade', 'subclade', 'qc.overallScore', 'qc.overallStatus', 'totalFrameShifts', 'frameShifts', 'aaInsertions')
  
  df3 <- df3[order(df3$Sample),]  # Ordenar por Sample (Hospital.sample.number)
  
  df3 <- df3[rowSums(is.na(df3)) < ncol(df3), ]
  
  write.table(df3, paste0(extra_dir, di, "_Nextclade_", RUN, ".csv"), row.names = F, sep=";")
  
  
  #--> !! En el report original puede haber clados asignados a muestras sin subtipo y dataset asignado, si tienen parte del HA.
  
  
  
  ### 4) Merge results ###
  dim(df1)
  dim(df2.cov)
  dim(df2.dep)
  dim(df2.total_map)
  dim(df3)
  
  # Combine map stats
  df2.def <- cbind(df2.cov, df2.dep[,2:9], df2.total_map$mapped_reads)
  colnames(df2.def)[18] <- "mapped_reads"
  
  # Merge tables by Sample name
  df4.0 <- merge(df1, df2.def, by.x="Sample", by.y="Sample", all = TRUE)
  df4 <- merge(df4.0, df3, by.x="Sample", by.y="Sample", all = TRUE)      # Merge with nextclade results (df3)
  
  
  ### 5) cols with constant values and empty cols ###
  nr <- nrow(df4)
  
  ID <- rep("", nr)
  SIP <- rep("", nr)
  Hospital <- rep("", nr)       # completar
  Reception.date <- rep("", nr) # completar
  Hospital.date <- rep("", nr)  # completar
  Relevance <- rep("VIGILANCIA", nr)
  Comments <- rep("", nr)
  Gender <- rep("", nr)
  Age <- rep("", nr)
  Status_Genomic <- rep("", nr) # Completar
  Status_HA.NA <- rep("", nr) # Completar
  runid <- rep(RUN, nr)
  
  trimmed_reads <- as.numeric(unlist(lapply(strsplit(df4$reads_after_trimming, " "), '[', 1))) # Formato correcto
  percent_missed_reads <- 100 - round(( as.numeric(df4$mapped_reads) / trimmed_reads)*100, 2)  # Calcular % missed reads
  
  version <- rep("Nextclade v3.13.0", nr)
  GISAID <- rep("", nr)
  
  ### 6) Ordenar columnas
  
  df5 <- cbind(ID, df4$Sample, SIP, Hospital, Reception.date, Hospital.date,
               Relevance, Comments, Gender, Age, Status_Genomic, Status_HA.NA, df4[,16:23],
               df4[24:31], runid, df4$IRMA_consensus_N_count, df4$reads_before_trimming,
               trimmed_reads, df4$mapped_reads, percent_missed_reads,
               df4[,5:7], df4[,8:11], df4$clade.y, df4$subclade.y, version, df4$qc.overallScore, 
               df4$qc.overallStatus, df4[,35:37], GISAID)
  
  colnames(df5) <- c("ID","Hospital.sample.number","SIP","Hospital","Reception.date","Hospital.date","Relevance","Comments","Gender","Age","Status_Genomic", "Status_HA.NA","HA.Coverage","MP.Coverage","NA.Coverage","NP.Coverage","NS.Coverage","PA.Coverage","PB1.Coverage","PB2.Coverage","HA.Mean_Depth","MP.Mean_Depth","NA.Mean_Depth","NP.Mean_Depth","NS.Mean_Depth","PA.Mean_Depth","PB1.Mean_Depth","PB2.Mean_Depth","runid","N_count","reads","reads_after_trimming","mapped_reads","percent_missed_reads","percent_Influenza.A","percent_Influenza.B","percent_Homo.sapiens","IRMA_type","abricate_InsaFlu_type","IRMA_subtype","abricate_InsaFlu_subtype","clade","subclade","version","qc.overallScore","qc.overallStatus","totalFrameShifts","frameShifts","aaInsertions","GISAID")
  
  
  ## Calcular Status de cada muestra
  
  # CRITERIOS POR DEFINIR
  #   * Failed si hay algún segmento con 0 cobertura
  #   * Failed si HA o NA tienen <90% cobertura
  
  #df5$Status[apply(is.na(df5[, 12:19]), FUN=any, MARGIN=1)] <- "Failed_QC" # Fail si al menos un segmento sin cobertura
  #df5$Status[is.na(df5$N_count) | df5$N_count > 3000] <- "Failed_QC"
  
  
  ## > QC status basado en cobertura sobre HA NA
  df5$HA.Coverage[is.na(df5$HA.Coverage)] <- 0
  df5$MP.Coverage[is.na(df5$MP.Coverage)] <- 0
  df5$NA.Coverage[is.na(df5$NA.Coverage)] <- 0
  df5$NP.Coverage[is.na(df5$NP.Coverage)] <- 0
  df5$NS.Coverage[is.na(df5$NS.Coverage)] <- 0
  df5$PA.Coverage[is.na(df5$PA.Coverage)] <- 0
  df5$PB1.Coverage[is.na(df5$PB1.Coverage)] <- 0
  df5$PB2.Coverage[is.na(df5$PB2.Coverage)] <- 0
  
  
  df5$Status_HA.NA[as.numeric(df5$HA.Coverage) < 90 & as.numeric(df5$NA.Coverage) < 90 ] <- "Failed"
  
  df5$Status_HA.NA[as.numeric(df5$HA.Coverage) < 90 & as.numeric(df5$NA.Coverage) >= 90 ] <- "Failed_HA"
  df5$Status_HA.NA[as.numeric(df5$HA.Coverage) >= 90 & as.numeric(df5$NA.Coverage) < 90 ] <- "Failed_NA"

  df5$Status_HA.NA[df5$Status_HA.NA == ""] <- "Sequenced"
  
  
  ## > QC status basado en cobertura del genoma completo.
  
  # Longitud de cada segmento:
  segment <- c("RSV_BD","RSV_AD","RSV_B","RSV_A","A_MP","A_NP","A_NS","A_PA","A_PB1","A_PB2","A_HA_H1","A_HA_H10","A_HA_H11","A_HA_H12","A_HA_H13","A_HA_H14","A_HA_H15","A_HA_H16","A_HA_H2","A_HA_H3","A_HA_H4","A_HA_H5","A_HA_H6","A_HA_H7","A_HA_H8","A_HA_H9","A_NA_N1","A_NA_N2","A_NA_N3","A_NA_N4","A_NA_N5","A_NA_N6","A_NA_N7","A_NA_N8","A_NA_N9","B_HA","B_MP","B_NA","B_NP","B_NS","B_PA","B_PB1","B_PB2")
  seg_len <- c(15277,15277,15225,15223,982,1497,863,2151,2274,2280,1704,1686,1698,1695,1701,1707,1713,1698,1689,1704,1695,1707,1704,1713,1701,1683,1413,1410,1410,1413,1422,1413,1416,1413,1413,1758,1139,1408,1683,1034,2181,2263,2313)
  
  
  names(seg_len) <- segment
  
  segment_A <- seg_len[5:10]
  segment_A["HA"] <- median(seg_len[11:26])
  segment_A["NA"] <- median(seg_len[27:35])
  
  names(segment_A) <- c("MP", "NP", "NS", "PA", "PB1", "PB2", "HA", "NA")
  
  segment_B <- seg_len[36:43]
  names(segment_B) <- segment_names
  
  
  
  ## 2) coverage, depth y mapped reads ###
  
  # indir <- "reports/"
  # afile <- "merged_bam_coverage_results.tsv"
  # 
  # f2 <- read.table(paste0(basedir, indir, afile), sep="\t", header = T, row.names = 1) # read table
  
  # !! --> No hay resultado de 1 negativo!!
  
  # Separar resultados coverage de mapped+depth
  #f2.mapped_dep <- f2[, 1:24]
  f2.cov_len   <- f2[, 25:dim(f2)[2]]
  f2.cov_ref   <- f2[, 25:dim(f2)[2]]
  
  # Get segment names
  segment_names <- c("HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2")
  df2 <- as.data.frame(rownames(f2))
  df2[,2:9] <- NA # create 8 empty columnes
  colnames(df2) <- c("Sample", segment_names)
  rownames(df2) <- df2$Sample
  
  # Tables for depth, coverage
  df2.cov_len <- df2
  df2.cov_ref <- df2

  
  ### Loop ids and segments <- get reference and assembly len
  for ( id in rownames(f2) ) {
    for (iseg in segment_names){
      
      # Get subtype:
      stype <- df5[df5$Hospital.sample.number==id, "IRMA_type"]
      
      
      # Coverage
      idx <- which(grepl(paste0("_",iseg), f2.cov_len[id,])) # grep idx 1 segment
      
      if (length(idx)==0){
        df2.cov_len[id, iseg] <- NA # Si no hay valores, introducir NA
      
      
        # Si no hay valores, introducir un tamaño de referencia medio. 
        if (stype == "Type_A") {
          df2.cov_ref[id, iseg] <- segment_A[iseg] # Si no hay valores, se introducen los de la referencia o una mediana
        } else if (stype == "Type_B") {
          df2.cov_ref[id, iseg] <- segment_B[iseg] # Si no hay valores, se introducen los de la referencia o una mediana
        } else {
          df2.cov_ref[id, iseg] <- segment_A[iseg] # Si no hay valores, se introducen los de la referencia o una mediana
        } 
      
      } else {
      value_co <- f2.cov_len[id, c(idx,(idx+2))]  # obtener los valores asociados a coverage
      df2.cov_len[id, iseg] <- value_co[2] # introducir valor percent_coverage
      
      value_co <- f2.cov_ref[id, c(idx,(idx+1))]  # obtener los valores asociados a coverage
      df2.cov_ref[id, iseg] <- value_co[2] # introducir valor percent_coverage
      
      } 
    }
  }

  
  # Tamaño completo genoma vs tamaño ensamblado; suma de filas (segmentos)
  df2.cov_len[is.na(df2.cov_len)] <- 0
  
  Assembly_len <- rowSums(df2.cov_len[2:9])
  Ref_len      <- rowSums(df2.cov_ref[2:9])
  
  Genomic_cov <- round((Assembly_len/Ref_len)*100, 2)
  
  tmp.GenomCov <- as.data.frame(cbind(Genomic_cov, names(Genomic_cov)))
  tmp.GenomCov <- tmp.GenomCov[order(tmp.GenomCov$V2),]
  colnames(tmp.GenomCov) <- c("Genomic_Cov", "Sample")
  
  # Merge genomic coverage with the table
  df5.0 <- merge(df5, tmp.GenomCov, by.x="Hospital.sample.number", by.y="Sample", all = TRUE)
  
  # 90% coverage check
  df5.0$Status_Genomic[as.numeric(df5.0$Genomic_Cov) < 90 ] <- "Failed"
  df5.0$Status_Genomic[is.na(df5.0$Genomic_Cov)] <- "Failed"
  #df5.0$Status_Genomic[df5.0$Status_Genomic == ""] <- "Sequenced"
  df5.0$Status_Genomic[df5.0$Genomic_Cov >= 90] <- "Sequenced"
  
  ## Save new report
  
  df6 <- cbind(ID, df5.0$Hospital.sample.number, SIP, Hospital, Reception.date, Hospital.date,
               Relevance, Comments, Gender, Age, df5.0$Status_Genomic, df5$Status_HA.NA, df5.0$Genomic_Cov, df4[,16:23],
               df4[24:31], runid, df4$IRMA_consensus_N_count, df4$reads_before_trimming,
               trimmed_reads, df4$mapped_reads, percent_missed_reads,
               df4[,5:7], df4[,8:11], df4$clade.y, df4$subclade.y, version, df4$qc.overallScore, 
               df4$qc.overallStatus, df4[,37:39], GISAID)
  
  colnames(df6) <- c("ID","Hospital.sample.number","SIP","Hospital","Reception.date","Hospital.date","Relevance","Comments","Gender","Age","Status_Genomic", "Status_HA.NA", "Genome.Coverage" ,"HA.Coverage","MP.Coverage","NA.Coverage","NP.Coverage","NS.Coverage","PA.Coverage","PB1.Coverage","PB2.Coverage","HA.Mean_Depth","MP.Mean_Depth","NA.Mean_Depth","NP.Mean_Depth","NS.Mean_Depth","PA.Mean_Depth","PB1.Mean_Depth","PB2.Mean_Depth","runid","N_count","reads","reads_after_trimming","mapped_reads","percent_missed_reads","percent_Influenza.A","percent_Influenza.B","percent_Homo.sapiens","IRMA_type","abricate_InsaFlu_type","IRMA_subtype","abricate_InsaFlu_subtype","clade","subclade","version","qc.overallScore","qc.overallStatus","totalFrameShifts","frameShifts","aaInsertions","GISAID")

  
  ### 8) Write table
  di <- format(Sys.time(), "%y%m%d") # date of analysis
  
  outfile = paste0(outdir0, di, "_RESULTS_", RUN, ".csv")
  write.table(df6, outfile, row.names = F, sep=";")
  
}
# 




