# VC Test for Microbiome Data

**VCmicrobiome.jl** is a Julia package for testing microbiome effect in longitudinal study. It utilizes three types of exact tests

* exact likelihood ratio test (eLRT)
* exact restricted likelihood ratio test (eRLRT)
* exact score test (eScore)

# What's new in VCmicrobiome.jl

* overall microbiome effect
* clustered microbiome effects
* longitudinal case
* parallel testing for multiple phenotypes

For testing the overall microbiome effect, the input files for _VCmicrobiome.jl_ are microbiome `kernel distance matrix (.csv file)`, `covariates file (.csv file)` and `phenotype file (.csv file)`. Since OTUs can be clustered at higher phylogenetic level, for testing multiple microbiome clusters, the kernel name list file is required.

You should have these input files prepared before running our program. The `output file (.out file)` is a simple comma-delimited file containing the _p_-values for each phenotype and each group of OTUs under a certain testing scheme (eLRT, eRLRT or eScore).

To use **VCmicrobiome.jl** , you need to call the `microvctest()` function.


## Contents

* [Installation](hhttp://vcmicrobiomejl.readthedocs.io/en/latest/Installation/)
* [Dataformats](http://vcmicrobiomejl.readthedocs.io/en/latest/Dataformats/)
* [Usage](http://vcmicrobiomejl.readthedocs.io/en/latest/Usage/)
* [Examples](http://vcmicrobiomejl.readthedocs.io/en/latest/Examples/)
