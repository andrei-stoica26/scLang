# scLang

`scLang` is a suite for package development for scRNA-seq analysis. 
It offers functions—both getters and setters—that can operate 
on both Seurat and SingleCellExperiment objects, intended to help 
developers build tools compatible with both types of input.

## Installation

To install `scLang`, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("scLang")
```
