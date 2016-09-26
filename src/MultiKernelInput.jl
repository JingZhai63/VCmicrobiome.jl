

function KernlInput(args...; KernelList::String = "",nObs::Int = 99) 
    Klist = readdlm(KernelList, ',')
    Klist = Klist[2:end,:]
    nKernel = length(Klist)
    function KernelArray(arg...)
	  KK=Array(Float64,nObs,nObs)
	  return KK
	end
    Ktuple=ntuple(length(Klist),KernelArray::Function);
    for i = 1:length(Klist)
    	 K= readdlm(Klist[i], ',')
    	 K=K[2:end,:]
         K=convert(Array{Float64, 2}, K)
         copy!(Ktuple[i],K)
    end
    return Ktuple
end

