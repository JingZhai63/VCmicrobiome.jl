# Installation

To install _VCmicrobiome.jl_, open up Julia and then type

```julia
julia> Pkg.clone("https://github.com/JingZhai63/VCmicrobiome.jl.git")
```
If the Julia packages _StatsBase_, _Rmath_ and _ConjugatePriors_ are not installed, please install them first

```julia
julia> Pkg.add("StatsBase")
julia> Pkg.add("Rmath")
julia> Pkg.build("Rmath")
import Rmath: libRmath
libRmath
julia> Pkg.add("ConjugatePriors")
```

Then you can use the following command to verify that the  _VCmicrobiome.jl_ package has been installed successfully

```julia
julia> using VCmicrobiome
julia> microvctest()
microvctest (generic function with 1 method)
```
