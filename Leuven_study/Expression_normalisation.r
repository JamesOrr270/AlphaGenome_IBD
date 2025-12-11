#This data has until the comparison of data. Loads the cell files, normalises using RMA and then annotates the probe IDs

library(oligo)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(dplyr)

setwd("/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/expression_data")
celFiles <- list.celfiles()
rawData <- read.celfiles(celFiles)

# RMA normalization
patient_manfiltered <- oligo::rma(rawData, target = 'core')

# Probe  annotation
anno_patient <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(patient_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_patient <- subset(anno_patient, !is.na(SYMBOL))

anno_grouped <- group_by(anno_patient, PROBEID)
anno_summarized <-
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)

head(anno_filtered)

probe_stats <- anno_filtered

nrow(probe_stats)

ids_to_exlude <- (featureNames(patient_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)

patient_final <- subset(patient_manfiltered, !ids_to_exlude)
validObject(patient_final)

head(anno_patient)

fData(patient_final)$PROBEID <- rownames(fData(patient_final))

fData(patient_final) <- left_join(fData(patient_final), anno_patient)

rownames(fData(patient_final)) <- fData(patient_final)$PROBEID
validObject(patient_final)

library(tibble)

expr_with_anno <- exprs(patient_final) %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  left_join(fData(patient_final), by = "PROBEID")

write.csv(expr_with_anno, 
          "patient_annotated_expression.csv", 
          row.names = FALSE)
