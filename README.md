# VC Test for Longitudinal Microbiome Study

**VCTestMicrobiome.jl** is a Julia package for performing exact variance component tests in longitudinal microbiome study. It utilizes three types of exact tests from Julia package **VarianceComponentTest.jl**

* exact likelihood ratio test (eLRT)
* exact restricted likelihood ratio test (eRLRT)
* exact score test (eScore)

The input files for _VCTestMicrobiome.jl_ are microbiome kernel distance matrix (.csv file), covariates file (.csv file) and phenotype file (.csv file). You should have these input files prepared before running our program. The output file (.out file) is a simple comma-delimited file containing the p-values for each phenotype and each group of OTUs under a certain testing scheme (eLRT, eRLRT or eScore).

To use VCTestMicrobiome.jl , you need to call the MicrobiomeVCTest() function.

Please see the detailed documents in [Read the Docs](http://127.0.0.1:8000) 

# Installation

To install _VCTestMicrobiome.jl_, open up Julia and then type

```julia
julia> Pkg.update()
julia> Pkg.clone("git://src/MicrobiomeVCTest.jl")
```
