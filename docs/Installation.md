# Installation

To install _VCTestMicrobiome.jl_, open up Julia and then type

```julia
julia> Pkg.clone("https://github.com/JingZhai63/VCmicrobiome.jl.git")
```
If the Julia package _VarianceComponentTest_ is not installed, please install first

```julia
julia> Pkg.add("VarianceComponentTest")
```

Then you can use the following command to verify that the package has been installed successfully

```julia
julia> using VCTestMicrobiome
julia> MicrobiomeVCTest()
MicrobiomeVCTest (generic function with 1 method)
```

Also, to install Tao-Hu's Julia package _VarianceComponentTest.jl_, please see [installation](http://variancecomponenttestjl.readthedocs.io/en/latest/installation/)
