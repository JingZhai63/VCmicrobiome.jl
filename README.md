# VC Test for Microbiome Data

VCTestMicrobiome.jl is a Julia package for performing exact variance component tests in longitudinal microbiome study. It provides three types of exact tests

* exact likelihood ratio test (eLRT)
* exact restricted likelihood ratio test (eRLRT)
* exact score test (eScore)

The input files for _VCTestMicrobiome.jl_ are microbiome kernel distance matrix (.csv file), covariates file (.csv file) and phenotype file (.csv file). You should have these input files prepared before running our program. The output file (.out file) is a simple comma-delimited file containing the p-values for each group of OTUs under a certain testing scheme (eLRT, eRLRT or eScore).

To use VCTestMicrobiome.jl , you need to call the MicrobiomeVCTest() function.


## Contents

* [Installation](http://127.0.0.1:8000/Installation/)
* [Dataformats](http://127.0.0.1:8000/Dataformats/)
* [Usage](http://127.0.0.1:8000/Usage/)
* [Examples](http://127.0.0.1:8000/Examples/)

