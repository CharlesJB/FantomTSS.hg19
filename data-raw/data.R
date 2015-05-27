# Prepare the datasets
#
# This function will create and save the datasets in the \code{R/systada.rda}
# file.
#
# param force If TRUE, the inst/extdata file will be created, even if it
#             exists. Default: \code{FALSE}.
#
# examples
#  prepare_internal_datasets
prepare_datasets <- function(force = FALSE) {
  # Add inst/extdata
  if (!dir.exists("inst/extdata/")) {
    dir.create("inst/extdata/", recursive = TRUE)
  }

  exp_description <<- data_fantom_exp_description(force = force)
  devtools::use_data(exp_description, internal = FALSE, overwrite = TRUE)
  metadata <- metadata_fantom_tss(force = force)
  get_sub_metadata <- function(i) {
    start <- (1+(i-1)*(ncol(metadata)/5))
    end <- (i*(ncol(metadata)/5))
    metadata[start:end]
  }
  metadata_1 <<- get_sub_metadata(1)
  metadata_2 <<- get_sub_metadata(2)
  metadata_3 <<- get_sub_metadata(3)
  metadata_4 <<- get_sub_metadata(4)
  metadata_5 <<- get_sub_metadata(5)
  devtools::use_data(metadata_1, metadata_2, metadata_3, metadata_4, metadata_5,
                    internal = FALSE, overwrite = TRUE)
  tss <<- data_fantom_tss(metadata = metadata)
  devtools::use_data(tss, internal = FALSE, overwrite = TRUE)
}

# Prepare the TSS metadata dataset
#
# This will download the Fantom's TSS file and convert it into \code{GRanges}
# format. The file is downloaded in the \code{inst/extdata/} directory.
#
# Download url:
#   "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/"
# Filename:
# "hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz"
#
# param force If TRUE, the inst/extdata file will be created, even if it
#             exists. Default: \code{FALSE}.
#
# return A \code{DataFrame} object with all the metadata
#
# examples
#   \dontrun{
#     FantomEnhancers.hg19::metadata_fantom_tss()
#   }
metadata_fantom_tss <- function(force = FALSE) {
  # Download TSS
  filename <- "hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz"
  url <- "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/"
  url <- paste0(url, filename)
  filename <- paste0("inst/extdata/", filename)
  download_file(url, filename, force = force)

  # Fetch header infos
  col_desc <- readLines(filename, 5000)
  col_desc <- col_desc[grepl("ColumnVariables", col_desc)]
  col_desc <-  gsub("\\[|\\]", ".", col_desc)
  col_desc <- strsplit(col_desc, "\\.")
  header <- sapply(col_desc, function(x) {
		   i <- grepl("CNhs", x)
		   ifelse(any(i), x[i][1], x[2])})
  # Read table
  col_desc <- readLines(filename, 5000)
  col_desc <- substring(col_desc, 1, 3)
  nskip <- min(which(col_desc == "chr")) - 1
  filename <- paste("zcat", filename)
  df <- fread(filename, sep = "\t", header = FALSE, skip = nskip)
  setnames(df, header)
  df <- S4Vectors::DataFrame(df)
  i <- grepl("CNhs", colnames(df))
  for (name in colnames(df)[i]) {
    df[[name]] <- S4Vectors::Rle(df[[name]])
  }
  df
}

# Prepare the TSS dataset
#
# param force If TRUE, the inst/extdata file will be created, even if it
#             exists. Default: \code{FALSE}.
# param metadata Metadatas produced with \code{metadata_fantom_tss}. If
#                \code{NULL}, will use call the function to generate it.
#
# return The \code{GRanges} produced.
#
# examples
#   \dontrun{
#     FantomEnhancers.hg19::data_fantom_tss()
#   }
data_fantom_tss <- function(metadata = NULL, force = FALSE) {
  if (is.null(metadata)) {
    metadata <- metadata_fantom_tss(force = force)
  }

  # Convert to GRanges
  gr <- data.frame(do.call("rbind", strsplit(metadata[,1], ":|\\.\\.|\\,")))
  setnames(gr, c("seqnames", "start", "end", "strand"))
  gr[["start"]] <- as.numeric(gr[["start"]])
  gr[["end"]] <- as.numeric(gr[["end"]])
  tmp <- gr[["start"]]
  i <- gr[["start"]] > gr[["end"]]
  gr[["start"]][i] <- gr[["end"]][i]
  gr[["end"]][i] <- tmp[i]
  seqinfo <- GenomeInfoDb::Seqinfo(genome = "hg19")
  GenomicRanges::makeGRangesFromDataFrame(gr, seqinfo = seqinfo)
}

# Prepare the exp_description dataset
#
# This will download the Fantom's experiment description file and import it as
# a data.frame. The file will be downloaded in the \code{inst/extdata/}
# directory.
#
# Download url:
#  "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/"
# Filename:
#  "inst/extdata/00_human.cell_line.hCAGE.hg19.assay_sdrf.txt"
#
# param force If TRUE, the inst/extdata file will be created, even if it
#             exists. Default: \code{FALSE}.
#
# return The \code{data.frame} produced.
#
# examples
# FantomEnhancers.hg19::data_fantom_enhancers()
data_fantom_exp_description <- function(force = FALSE) {
  # Download TSS
  filename <- "00_human.cell_line.hCAGE.hg19.assay_sdrf.txt"
  url <-
    "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/"
  url <- paste0(url, filename)
  filename <- paste0("inst/extdata/", filename)
  download_file(url, filename, force = force)

  # Prepare the data.frame
  read.table(filename, header = TRUE, sep = "\t", quote = "",
             stringsAsFactors = FALSE)
}
