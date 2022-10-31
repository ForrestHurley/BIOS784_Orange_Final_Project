## Download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99254
## Paper: https://doi.org/10.1038/s41591-018-0045-3
library(SummarizedExperiment)
## Download data
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99254/suppl/GSE99254_NSCLC.TCell.S12346.count.txt.gz",
  destfile = "data/source_data/GSE99254_NSCLC.TCell.S12346.count.txt.gz",
  skip = T
)
R.utils::gunzip("data/source_data/GSE99254_NSCLC.TCell.S12346.count.txt.gz",
                skip = T)

## Download sample info
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99254/soft/GSE99254_family.soft.gz",
  destfile = "data/source_data/GSE99254_family.soft.gz",
  skip = T
)
R.utils::gunzip("data/source_data/GSE99254_family.soft.gz",
                skip = T)


## Read data and create summarized experiment
d <- data.table::fread("data/source_data/GSE99254_NSCLC.TCell.S12346.count.txt")
rownames(d) <- d$geneID

info <- data.table::fread(cmd = "grep -E -v '!|\\^|#' data/source_data/GSE99254_family.soft")

se <-
  SummarizedExperiment(assays = list(counts = d[, info$UniqueCell_ID, with = F]),
                       rowData = d[, c(1, 2)],
                       colData = info)

saveRDS(se, "data/processed_data/se.rds")