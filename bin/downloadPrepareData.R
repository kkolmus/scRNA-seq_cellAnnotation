### LIB ###

suppressPackageStartupMessages({
  library(futile.logger)
  library(GEOquery)
  library(R.utils)
  library(dplyr)
  library(tibble)
  library(stringr)
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
file.remove(file.path(raw.data.dir, GEOid, "GSE127465_human_counts_normalized_54773x41861.mtx.gz"), 
            file.path(raw.data.dir, GEOid, "GSE127465_mouse_counts_normalized_15939x28205.mtx.gz"))


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


ls.metadata <- list.files(path = file.path(raw.data.dir, GEOid), pattern = "GSE")
ls.metadata <- ls.metadata[-5]

for (f in ls.metadata) {
  flog.debug(paste0("Unpacking the following file: ", f))
  gunzip(filename = file.path(raw.data.dir, GEOid, f))
}

ls.metadata <- list.files(path = file.path(raw.data.dir, GEOid), pattern = "GSE")
ls.metadata <- ls.metadata[-5]

for (f in ls.metadata) {
  print(f)
  temp.file.name <- substring(f, 11, (nchar(f)-4))
  temp.file <- read.table(file.path(raw.data.dir, GEOid, f), fill = TRUE, header = TRUE, sep = "\t")
  assign(temp.file.name, temp.file)
}

saveRDS(human_cell_metadata_54773x25, file.path(data.input.dir, "human.metadata.RDS"))
saveRDS(mouse_cell_metadata_15939x12, file.path(data.input.dir, "mouse.metadata.RDS"))

# subset count matrix to those cells that are annotated by authors


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

samples.ls <- list.files(
  path = file.path(raw.data.dir, GEOid),
  pattern = "_raw_counts.tsv")


flog.debug("Count data for human samples")

human.samples.ls <- samples.ls[str_detect(samples.ls, pattern = "human")]
human.samples.id <- c()
for (f in human.samples.ls) {
  temp.file.name = substring(f, 1, 10)
  human.samples.id <- c(human.samples.id, temp.file.name)
}

human.samples.id <- setNames(vector("list", length(human.samples.id)), human.samples.id)

for (f in human.samples.ls) {
  temp.file.name = substring(f, 1, 10)
  flog.debug(paste("Processing sample", temp.file.name, sep = " "))
  temp.file <- read.table(file.path(raw.data.dir, GEOid, f))
  temp.file <- t(temp.file)
  temp.file <- header.true(temp.file)
  rownames(temp.file) <- NULL
  temp.file <- column_to_rownames(temp.file, "barcode")
  temp.colnames <- paste0(temp.file.name, "_", colnames(temp.file))
  colnames(temp.file) <- temp.colnames
  temp.file <- as.matrix(temp.file)
  saveRDS(temp.file, file.path(data.input.dir, paste0(temp.file.name, ".RDS")))
  human.samples.id[[temp.file.name]] <- temp.file
  rm(temp.file)
}

saveRDS(human.samples.id, file.path(data.input.dir, "ls.countMatrix.human.RDS"))


flog.debug("Master count matrix for human samples")

df <- human.samples.id[[1]]
for (i in 2:length(human.samples.id)) {
  df <- merge(df, human.samples.id[[i]], by = row.names, all = TRUE)
  rownames(df) <- df$Row.names
  df <- df[, !(names(df) %in% "Row.names")]
}

human.count.matrix <- as.matrix(df)

saveRDS(human.count.matrix, file.path(data.input.dir, "countMatrix.human.RDS"))


flog.debug("Count data for mouse samples")

mouse.samples.ls <- samples.ls[str_detect(samples.ls, pattern = "mouse")]
mouse.samples.id <- c()
for (f in mouse.samples.ls) {
  temp.file.name = substring(f, 1, 10)
  mouse.samples.id <- c(mouse.samples.id, temp.file.name)
}

mouse.samples.id <- setNames(vector("list", length(mouse.samples.id)), mouse.samples.id)

for (f in mouse.samples.id) {
  temp.file.name = substring(f, 1, 10)
  flog.debug(paste("Processing sample", temp.file.name, sep = " "))
  temp.file <- read.table(file.path(raw.data.dir, GEOid, f))
  temp.file <- t(temp.file)
  temp.file <- header.true(temp.file)
  rownames(temp.file) <- NULL
  temp.file <- column_to_rownames(temp.file, "barcode")
  temp.colnames <- paste0(temp.file.name, "_", colnames(temp.file))
  colnames(temp.file) <- temp.colnames
  temp.file <- as.matrix(temp.file)
  saveRDS(temp.file, file.path(data.input.dir, paste0(temp.file.name, ".RDS")))
  mouse.samples.id[[temp.file.name]] <- temp.file
  rm(temp.file)
}

saveRDS(mouse.samples.id, file.path(data.input.dir, "ls.countMatrix.mouse.RDS"))


flog.debug("Master count matrix for mouse samples")

df <- mouse.samples.id[[1]]
for (i in 2:length(mouse.samples.id)) {
  df <- merge(df, mouse.samples.id[[i]], by = row.names, all = TRUE)
  rownames(df) <- df$Row.names
  df <- df[, !(names(df) %in% "Row.names")]
}

mouse.count.matrix <- as.matrix(df)

saveRDS(mouse.count.matrix, file.path(data.input.dir, "countMatrix.mouse.RDS"))


### SESSION INFO ###
flog.debug("SessionInfo")
sessionInfo()