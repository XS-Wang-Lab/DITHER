# DITHER
An algorithm for Defining IntraTumor Heterogeneity based on EntRopy.
## Description
The package is used to calculate IntraTumor Heterogeneity (ITH) based on somatic mutation and copy number alteration (CNA) profiles in tumors.
## Details
- The input of the function ```mDITHER()``` should be a tab-delimited plain text file, which should include at least two columns ("Sample" and "VAF"). The "Sample" column is the tumor sample ID; the "VAF" column is the somatically mutant allele fractions (MAFs) among loci in the tumor.
- The input of the function ```cDITHER()``` should be a tab-delimited plain text file, which should include at least two columns ("Sample" and "Segment_Mean"). The "Sample" column is the tumor sample ID; the "Segment_Mean" column is the mean segments in the tumor, which equal to log2(copy number/2).
- The input of the function ```DITHER()``` should include the two files mentioned above.
## Installation
You can install the released version of **DITHER** with:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("XS-Wang-Lab/DITHER/DITHER")
```
## Examples
```r
## Calculate ITH based on the entropy of somatic mutation profiles in tumors ----------
library(DITHER)
path = system.file("extdata", "example_mDITHER.txt", package = "DITHER", mustWork = TRUE)
input_data_mut = read.table(path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
mDITHER(input_data_mut)
```

```r
## Calculate ITH based on the entropy of CNA profiles in tumors ----------
library(DITHER)
path = system.file("extdata", "example_cDITHER.txt", package = "DITHER", mustWork = TRUE)
input_data_cna = read.table(path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
cDITHER(input_data_cna)
```

```r
## Calculate ITH based on the entropy of somatic mutation and CNA profiles in tumors ----------
library(DITHER)
path1 = system.file("extdata", "example_mDITHER.txt", package = "DITHER", mustWork = TRUE)
path2 = system.file("extdata", "example_cDITHER.txt", package = "DITHER", mustWork = TRUE)
input_data_mut = read.table(path1, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
input_data_cna = read.table(path2, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
DITHER(input_data_mut, input_data_cna)
```

## Citation
Li L, Chen C, Wang X. DITHER: an algorithm for Defining IntraTumor Heterogeneity based on EntRopy. Brief Bioinform. 2021 Nov 5;22(6):bbab202. doi: 10.1093/bib/bbab202. PMID: 34096997.
