__precompile__()

module VCmicrobiome


using ConjugatePriors
using Rmath
using StatsBase
using Distributions

export

      vctest,
      vctestnullsim,
      microvctest

include("vctest.jl")
include("vctestnullsim.jl")
include("microvctest.jl")
# include("MultiKernelInput.jl")
# include("MMicrobiomeVCTest.jl")
end # module
