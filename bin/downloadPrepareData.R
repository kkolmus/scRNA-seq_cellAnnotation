### LIB ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(GEOquery)
  library(R.utils)
  library(dplyr)
  library(tibble)
})

### FUNCTION ###
header.true <- function(input) {
  df <- as.data.frame(input)
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}


### MAIN ###

flog.threshold(DEBUG)

flog.debug("Download scRNAseq data")

GEOid <- "GSE127465"
raw.data.dir <- "~/Desktop/HT projects/scRNA-seq_cellAnnotation/data/raw.data"
dir.create(raw.data.dir)
data.input.dir <- "~/Desktop/HT projects/scRNA-seq_cellAnnotation/data/input"
dir.create(data.input.dir)
data.output.dir <- "~/Desktop/HT projects/scRNA-seq_cellAnnotation/data/output"
dir.create(data.output.dir)

gds <- getGEO(GEO = GEOid, destdir = raw.data.dir, GSEMatrix = TRUE)
getGEOSuppFiles(GEO = GEOid, baseDir = raw.data.dir) # supplementray file download
# getGEOfile(GEO = GEOid, destdir = raw.data.dir) # soft format files download


flog.debug("Prepare data input for downstream analysis")

flog.debug("Metadata")
pdata <- lapply(gds, phenoData)
pdata.human <- pdata[["GSE127465-GPL18573_series_matrix.txt.gz"]]@data
pdata.human <- pdata.human[, c(2,9,8,52:62)]
head(pdata.human)
colnames(pdata.human) <- c("geo_accession", "organism", "source", "age", "gender", "cancer_type",
                           "PD-L1 expression", "purification_method", "mortality", "tissue_diagnosis",
                           "treatment_prior_to_surgery", "tumor_stage", "surger_type", "treatment_type")
saveRDS(pdata.human, file.path(data.input.dir, "pdata.human.RDS"))
pdata.mouse <- pdata[["GSE127465-GPL19057_series_matrix.txt.gz"]]@data
pdata.mouse <- pdata.mouse[, c(2,9,55,8,57,54,51,52,53)]
colnames(pdata.mouse) <- c("geo_accession", "organism", "strain", "source", "tumor_cell_line",  
                           "purification", "age_at_tumor_injection", "age_at_sacrifice", "gender")
saveRDS(pdata.mouse, file.path(data.input.dir, "pdata.mouse.RDS"))

flog.debug("Count data")
# untar data files

# system("defaults write org.R-project.R force.LANG en_US.UTF-8")
# untar(tarfile = file.path(raw.data.dir, GEOid, "GSE127465_RAW.tar"), list = TRUE)
untar(tarfile = file.path(raw.data.dir, GEOid, "GSE127465_RAW.tar"),
      exdir = raw.data.dir)
# unzip data files
count.data.files <- list.files(file.path(raw.data.dir, GEOid), 
                               pattern = "raw_counts.tsv.gz", 
                               full.names = FALSE)
for (f in count.data.files) {
  flog.debug(paste0("Unpacking the following file: ", f))
  gunzip(filename = file.path(raw.data.dir, GEOid, f))
}

# test
data <- read.table(file.path(raw.data.dir, GEOid, "GSM3635278_human_p1t1_raw_counts.tsv"))
data <- t(data)
data <- header.true(data)
rownames(data) <- NULL
data <- column_to_rownames(data, "barcode")
data2 <- data[, 1:25]
# filter count data frames to select those cells with enough features
# understand files with metadata, do they how annotated cells?

### SESSION INFO ###
sessionInfo()