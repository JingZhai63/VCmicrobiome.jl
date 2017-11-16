module VCmicrobiome

import Distributions, VarianceComponentTest
using Distributions
using VarianceComponentTest
export

      MicrobiomeVCTest,
      KernlInput,
      MMicrobiomeVCTest

include("MicrobiomeVCTest.jl")
include("MultiKernelInput.jl")
# include("readNewick.jl")
include("MMicrobiomeVCTest.jl")
end # module
