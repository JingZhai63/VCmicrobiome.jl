# VC Test for Microbiome Data

**VCTestMicrobiome.jl** is a Julia package for testing microbiome effect in longitudinal study. It utilizes three types of exact tests from Tao-Hu's Julia package [VarianceComponentTest.jl](https://github.com/Tao-Hu/VarianceComponentTest.jl/blob/master/README.md)

* exact likelihood ratio test (eLRT)
* exact restricted likelihood ratio test (eRLRT)
* exact score test (eScore)

#What's new in VCTestMicrobiome.jl

* overall microbiome effect
* clustered microbiome effects
* longitudinal case
* parallel testing for multiple phenotypes

For testing the overall microbiome effect, the input files for _VCTestMicrobiome.jl_ are microbiome `kernel distance matrix (.csv file)`, `covariates file (.csv file)` and `phenotype file (.csv file)`. Since OTUs can be clustered at higher phylogenetic level, for testing multiple microbiome clusters, the kernel name list file is required. 

You should have these input files prepared before running our program. The `output file (.out file)` is a simple comma-delimited file containing the _p_-values for each phenotype and each group of OTUs under a certain testing scheme (eLRT, eRLRT or eScore).

To use **VCTestMicrobiome.jl** , you need to call the `MicrobiomeVCTest()` function.


## Contents

* [Installation](http://127.0.0.1:8000/Installation/)
* [Dataformats](http://127.0.0.1:8000/Dataformats/)
* [Usage](http://127.0.0.1:8000/Usage/)
* [Examples](http://127.0.0.1:8000/Examples/)

