# FantomTSS.hg19

The goal of this package is to easily extract the genomic positions and the normalized expression of Fantom's TSS (in TPM) in hg19.

## Installation

```
require(devtools)
devtools::install_github("CharlesJB/FantomTSS.hg19")
```

## Main functions

The `get_fantom_tss` functions returns a `GRanges` object with all the tss and no metadata columns:
```
get_fantom_tss()
```

The `get_fantom_tss_tpm` returns a `GRanges` object with all the tss and selected metadata columns:
```
# To get the expression of tss in A549
get_fantom_tss_tpm(cell_lines = "A549")

# To get the expression of tss in A549 and K562
get_fantom_tss_tpm(cell_lines = c("A549", "K562"))

# To merge all the columns for a cell line by calculating the mean value for
# each tss
get_fantom_tss_tpm(cell_lines = c("A549", "K562"), merge.FUN = mean)

# To merge all the columns for a cell line by calculating the sum for each
# tss
get_fantom_tss_tpm(cell_lines = c("A549", "K562"), merge.FUN = sum)
```

## Reproduce the datasets

If you wish to reproduce the datasets, you can use the `prepare_datasets` function. From a R session in the packages main directory:

```
library(devtools)
library(data.table)
load_all()
source("data-raw/data.R")
prepare_datasets()
```

Note that this will download large files (~800 MB) from Fantom's web page. You will also need to reinstall the package if you want to use those new datasets.
