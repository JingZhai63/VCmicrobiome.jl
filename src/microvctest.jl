using ConjugatePriors
using Rmath
using StatsBase

function microvctest(arg...; covFile::String = "", responseFile::String = "",
                     kernelFile::String = "",  kadjFile::String = "",
                     outFile::String = "", zFile::String = "", yIdx::Int64 = 3,
                     IDidx::Int64 = 1,  xIdx::Array{Int64, 1} = vcat(3,4,5),
                     out::Bool = false,
                     longitudinal::Bool = true, fine::Bool = false,
                     test::String = "eRLRT", device::String = "CPU",
                     infLambda::Float64 = 0.0, nMMmax::Int64 = 0,
                     nBlockAscent::Int64 = 1000, nNullSimPts::Int64 = 1000000,
                     nNullSimNewtonIter::Int64 = 15, tolX::Float64 = 1e-4,
                     MemoryLimit::Int64 = 200000000, lowRank::Float64 = 0.4,
                     pvalueComputing::String = "chi2")

##VCmicrobiome::--------------set default maximum MM iteration-----------------#
if nMMmax == 0  && test == "eLRT"
  nMMmax = 10
elseif nMMmax == 0 && test == "eRLRT"
  nMMmax = 1000
end

##VCmicrobiome::-----------------------get kernel------------------------------#
if isempty(kernelFile)
  error("mVctest:kernelwrongn\n", "# no kernel matrix")
else
  K = readdlm(kernelFile,',')
  if typeof(K[1,1]) != Float64
    K = K[2:end,:]
    K = convert(Array{Float64, 2}, K)
  else
    K = convert(Array{Float64, 2}, K)
  end
end

nObs = length(K[:, 1])

##VCmicrobiome::-----------------------get kernel adj -------------------------#
if fine & isempty(kadjFile)
  error("mVctest:kernelwrongn\n", "# no kernel adjustment matrix")
elseif fine & !isempty(kadjFile)
  Kadj = readdlm(kadjFile,',')
  if typeof(Kadj[1,1]) != Float64
    Kadj = Kadj[2:end,:]
    Kadj = convert(Array{Float64, 2}, Kadj)
  else
    Kadj = convert(Array{Float64, 2}, Kadj)
  end
end

#VCmicrobiome::------------------ get covariate and response-------------------#
if isempty(covFile)
  X = ones(nObs, 1) # N-by-1 matrix with all elements equal to 1
  X = convert(Array{Float64, 2}, X)
else
  covariate = readdlm(covFile,',')
  X = covariate[2:end, xIdx]
  if size(X, 1) != nObs
    error("VCmicrobiome:covwrongn\n",
    "# N in covariate file does not match dimension of kernel matrix")
  end
  X = convert(Array{Float64, 2}, X)
end

if isempty(responseFile)
  error("VCmicrobiome:noresponse\n", "# need to provide response file")
else
  Y = readdlm(responseFile,',')
  yname = Y[1, yIdx]
  Y = Y[2:end, yIdx]
  if size(Y, 1) != nObs
    error("VCmicrobiome:ywrongn\n",
    "# N in response file does not match dimension of Kernel")
  end
  y = convert(Array{Float64, 1}, Y[:])
end

if longitudinal && isempty(zFile)
  if isempty(IDidx) || isempty(covFile)
    error("VCmicrobiome:noid\n",
    "# need subject ID if Z matrix not provided")
  else
    id = covariate[2:end, IDidx]
  end
  dic = StatsBase.countmap(id)
  IDunq = unique(dic)
  Z = cat([1, 2],ones(dic[id[1]], 1))
  for i = 2:length(IDunq)
    Z = cat([1, 2], Z, ones(dic[id[i]], 1))
  end
  ZZ = Z * Z'
end

if isempty(zFile) == false
    Z = readdlm(zFile,',')
    Z = Z[2:end, :]
    Z = convert(Array{Float64, 2}, Z)
    ZZ = Z * Z'
end

offset = 0

# prepare output
if isempty(outFile)
  outFile = string(responseFile, "-", test,"-julia.out");
end

outresults = Array{Float64}(1, 3)
Xsvd = svdfact(X, thin = false);
rankX = countnz(Xsvd[:S] .> nObs * eps(Xsvd[:S][1]));
XtNullBasis = Xsvd[:U][:, rankX + 1 : end];

if longitudinal == false && fine == false
  ynew = zeros(nObs - rankX);
  BLAS.gemv!('N', 1.0, XtNullBasis', y, 0.0, ynew);
  yShift = Float64[];
  PrePartialSumW = [Float64[] Float64[]];
  PreTotalSumW = [Float64[] Float64[]];
  tmpn = nObs - rankX; ##zj::
  nPreRank = tmpn - 1
  partialSumWConst = Array{Float64}(nNullSimPts);
  totalSumWConst = Array{Float64}(nNullSimPts);
elseif longitudinal && fine == false
  QX = Xsvd[:U][:, 1 : rankX]
  QZsvd = svdfact(BLAS.gemm('T', 'N', XtNullBasis, Z))

  ZKeepIdx = (1 - QZsvd[:S]) .< 1e-6;
  rankQZ = countnz(ZKeepIdx)
  QZ = Array{Float64}(nObs, rankQZ);
  #orthonormal basis Q1 of C(Phi)-C(X)
  BLAS.gemm!('N', 'N', 1.0, Z, QZsvd[:V][:, 1 : rankQZ], 0.0, QZ);
  XPhitNullBasis = nullspace([QX Z]');
  tmpMat = Array{Float64}(rankQZ, size(Z,2));
  BLAS.gemm!('T', 'N', 1.0, QZ, Z, 0.0, tmpMat);
  PhiAdjsvd = svdfact(tmpMat, thin = false); #PhiAdjsvd ---WLambdaW'
  #W = PhiAdjsvd[:U]
  evalPhiAdj = PhiAdjsvd[:S] .^ 2
  # enforce first entry of each eigenvector to be >=0
  idxW = vec(PhiAdjsvd[:U][1, :]) .< 0
  PhiAdjsvd[:U][:, idxW] = - PhiAdjsvd[:U][:, idxW]

  KPhiAdj = PhiAdjsvd[:U][:, :]
  scale!(KPhiAdj, sqrt(evalPhiAdj / minimum(evalPhiAdj) - 1))
  InvSqrtevalPhiAdj = similar(evalPhiAdj)
  for i = 1 : length(evalPhiAdj)
    InvSqrtevalPhiAdj[i] = 1 / sqrt(evalPhiAdj[i])
  end
  weightedW = PhiAdjsvd[:U][:, :]
  scale!(weightedW, InvSqrtevalPhiAdj)  ##weightedW---W'*Lambda^1/2
  yShift = QZ' * y   # precompute shift in Y

  # prepare for simulating Chi Squares and sums of Chi Squares
  tmpn = length(evalPhiAdj)
  nPreRank = tmpn -1
  QRes = Array{Float64}(nObs, rankQZ) # Q2
  tmpvec = similar(yShift)
  tmpvecQRes = Array{Float64}(rankQZ)
  yWork = similar(evalPhiAdj)
  partialSumWConst = Array{Float64}(nNullSimPts)
  totalSumWConst = Array{Float64}(nNullSimPts)
  subXPhitSV = Array{Float64}(size(XPhitNullBasis, 2), rankQZ)
  pSubXPhitSV = pointer(subXPhitSV)
  tmpMat = Array{Float64}(nObs, size(XPhitNullBasis, 2))
  VWorkSqrt = Array{Float64}(length(evalPhiAdj), nObs)
  VWorkSqrt2 = Array{Float64}(length(evalPhiAdj), nObs)


  BLAS.gemm!('T', 'N', 1.0, K, XPhitNullBasis, 0.0, tmpMat)
  (UXPhitS, svalXPhitS, VXPhitS) = svd(tmpMat, thin = false)
  rankXPhitS = countnz(svalXPhitS .> size(XPhitNullBasis, 2) * eps(svalXPhitS[1]))
  pXPhitSV = pointer(VXPhitS)
  BLAS.blascopy!(size(XPhitNullBasis, 2) * rankQZ, pXPhitSV, 1, pSubXPhitSV, 1)
  BLAS.gemm!('N', 'N', 1.0, XPhitNullBasis, subXPhitSV, 0.0, QRes) # QRes---Q2
  #working Y and working non-trivial variance component
  BLAS.blascopy!(length(yShift), yShift, 1, tmpvec, 1) # tmpvec---Q1'*y
  BLAS.gemv!('T', 1.0, QRes, y, 0.0, tmpvecQRes)  # tmpvecQRes---Q2'*y
  BLAS.gemv!('N', 1.0, KPhiAdj, tmpvecQRes, 1.0, tmpvec) #tmpvec---Q1'*y+K*Q2'*Y
  BLAS.gemv!('T', 1.0, weightedW, tmpvec, 0.0, yWork)
  BLAS.gemm!('T', 'N', 1.0, QZ, K, 0.0, VWorkSqrt2)  #VWorkSqrt2=Q1'*Kernel
  BLAS.gemm!('T', 'N', 1.0, weightedW, VWorkSqrt2, 0.0, VWorkSqrt)
  Vform = "half"
elseif fine
  QX = Xsvd[:U][:, 1 : rankX];  #99 4 #orthonormal basis of Q0 of C(X)
  Phieig = eigfact(Kadj);

  if typeof(Phieig[:values][1]) == Complex{Float64}
    Phieigval = reinterpret(Float64, Phieig[:values][1:2:end])
  else
    Phieigval = reinterpret(Float64, Phieig[:values][1:end])
  end

  Phieigvec = zeros(size(Phieig[:vectors]))
  if typeof(Phieig[:vectors][1,1]) == Complex{Float64}
    for p in 1:size(Phieigvec,1)
        Phieigvec[p,:] = reinterpret(Float64, Phieig[:vectors][p,:][1:2:end])
    end
  else
    for p in 1:size(Phieigvec,1)
        Phieigvec[p,:] = reinterpret(Float64, Phieig[:vectors][p,:][1:end])
    end
  end

  rankPhi = countnz(Phieigval .> max(0.0, nPerKeep * eps(maximum(Phieigval))));
  idxEvecPhi = vec(Phieigvec[1, :]) .< 0;
  scale!(Phieigvec[:, idxEvecPhi], -1.0);
  rankPhi = min(rankPhi, floor(rankPhi * lowRank));
  # low rank approximation
  sortIdx = sortperm(Phieigval, rev = true);
  evalPhi = Phieigval[sortIdx[1 : Int(rankPhi)]];
  evecPhi = Phieigvec[:, sortIdx[1 : Int(rankPhi)]];
  QPhisvd = svdfact(BLAS.gemm('T', 'N', XtNullBasis, evecPhi));
  PhiKeepIdx = (1 - QPhisvd[:S]) .< 1e-6;
  rankQPhi = countnz(PhiKeepIdx);
  #QPhi = evecPhi * QPhisvd[:V][:, 1 : rankQPhi];
  QPhi = Array{Float64}(nPerKeep, rankQPhi);
  BLAS.gemm!('N', 'N', 1.0, evecPhi, QPhisvd[:V][:, 1 : rankQPhi], 0.0, QPhi);
  XPhitNullBasis = nullspace([QX evecPhi]');
  tmpMat = Array{Float64}(rankQPhi, Int(rankPhi));

  BLAS.gemm!('T', 'N', 1.0, QPhi, evecPhi, 0.0, tmpMat);
  scale!(tmpMat, sqrt(evalPhi));
  PhiAdjsvd = svdfact(tmpMat, thin = false);
  evalPhiAdj = PhiAdjsvd[:S] .^ 2;
  idxW = vec(PhiAdjsvd[:U][1, :]) .< 0;
  scale!(PhiAdjsvd[:U][:, idxW], -1.0);
  KPhiAdj = PhiAdjsvd[:U][:, :];
  scale!(KPhiAdj, sqrt(evalPhiAdj / minimum(evalPhiAdj) - 1));
  InvSqrtevalPhiAdj = similar(evalPhiAdj);
  for i = 1 : length(evalPhiAdj)
    InvSqrtevalPhiAdj[i] = 1 / sqrt(evalPhiAdj[i]);
  end
  weightedW = PhiAdjsvd[:U][:, :];
  scale!(weightedW, InvSqrtevalPhiAdj);
  yShift = QPhi' * y;

  # prepare for simulating Chi Squares and sums of Chi Squares
  tmpn = length(evalPhiAdj);
  nPreRank = length(evalPhiAdj)-1
  QRes = Array{Float64}(nPerKeep, rankQPhi);

  tmpvec = similar(yShift);
  tmpvecQRes = Array{Float64}(rankQPhi);
  yWork = similar(evalPhiAdj);
  partialSumWConst = Array{Float64}(nNullSimPts);
  totalSumWConst = Array{Float64}(nNullSimPts);
  subXPhitSV = Array{Float64}(size(XPhitNullBasis, 2), rankQPhi);
  pSubXPhitSV = pointer(subXPhitSV);

  tmpMat = Array{Float64}(nPerKeep, size(XPhitNullBasis, 2));
  VWorkSqrt = Array{Float64}(length(evalPhiAdj), nPerKeep);
  VWorkSqrt2 = Array{Float64}(length(evalPhiAdj), nPerKeep);

     # obtain some orthonormal vectors of space R^n - [QX,QPHI,S]
     # See Golub and Van Loan (1996) Algorithm 12.4.2 on p602
  BLAS.gemm!('T', 'N', 1.0, K, XPhitNullBasis, 0.0, tmpMat);
  (UXPhitS, svalXPhitS, VXPhitS) = svd(tmpMat, thin = false);

  pXPhitSV = pointer(VXPhitS)
  BLAS.blascopy!(size(XPhitNullBasis, 2) * rankQPhi, pXPhitSV, 1, pSubXPhitSV, 1);
  BLAS.gemm!('N', 'N', 1.0, XPhitNullBasis, subXPhitSV, 0.0, QRes); # QRes---Q2
  BLAS.blascopy!(length(yShift), yShift, 1, tmpvec, 1); # tmpvec---Q1'*y
  BLAS.gemv!('T', 1.0, QRes, y, 0.0, tmpvecQRes);  # tmpvecQRes---Q2'*y
  BLAS.gemv!('N', 1.0, KPhiAdj, tmpvecQRes, 1.0, tmpvec); #tmpvec---Q1'*y+K*Q2'*Y
  BLAS.gemv!('T', 1.0, weightedW, tmpvec, 0.0, yWork);

  BLAS.gemm!('T', 'N', 1.0, QPhi, K, 0.0, VWorkSqrt2);  #VWorkSqrt2=Q1'*Kernel
  BLAS.gemm!('T', 'N', 1.0, weightedW, VWorkSqrt2, 0.0, VWorkSqrt);
  Vform="half"
end

if pvalueComputing == "chi2"
  nPtsChi2 = 300;
  partialSumW = Array{Float64}(nPtsChi2);
  totalSumW = Array{Float64}(nPtsChi2);
  lambda = Array{Float64}(1, nPtsChi2);
  W = Array{Float64}(tmpn, nPtsChi2);
  tmpmat0 = Array{Float64}(tmpn, nPtsChi2);
  tmpmat1 = Array{Float64}(tmpn, nPtsChi2);
  tmpmat2 = Array{Float64}(tmpn, nPtsChi2);
  tmpmat3 = Array{Float64}(tmpn, nPtsChi2);
  tmpmat4 = Array{Float64}(tmpn, nPtsChi2);
  tmpmat5 = Array{Float64}(tmpn, nPtsChi2);
  denomvec = Array{Float64}(nPtsChi2);
  d1f = Array{Float64}(nPtsChi2);
  d2f = Array{Float64}(nPtsChi2);
  if test == "eScore"
    simnull = Array{Float64}(nPtsChi2);
  else
    simnull = Float64[];
  end
else
  nPtsChi2 = 300;
  partialSumW = Array{Float64}(nNullSimPts);
  totalSumW = Array{Float64}(nNullSimPts);
  lambda = Array{Float64}(1, nNullSimPts);
  W = Array{Float64}(tmpn, nNullSimPts);
  tmpmat0 = Array{Float64}(tmpn, nNullSimPts);
  tmpmat1 = Array{Float64}(tmpn, nNullSimPts);
  tmpmat2 = Array{Float64}(tmpn, nNullSimPts);
  tmpmat3 = Array{Float64}(tmpn, nNullSimPts);
  tmpmat4 = Array{Float64}(tmpn, nNullSimPts);
  tmpmat5 = Array{Float64}(tmpn, nNullSimPts);
  denomvec = Array{Float64}(nNullSimPts);
  d1f = Array{Float64}(nNullSimPts);
  d2f = Array{Float64}(nNullSimPts);
  if test == "eScore"
    simnull = Array{Float64}(nNullSimPts);
  else
    simnull = Float64[];
  end
end

WPreSim = randn(nObs, nNullSimPts);
for idxWc = 1 : nNullSimPts
  for idxWr = 1 : nObs
    WPreSim[idxWr, idxWc] = WPreSim[idxWr, idxWc] * WPreSim[idxWr, idxWc];
  end
end

  # simulate Chi Squares and sums of Chi Squares
PrePartialSumW = Array{Float64}(nNullSimPts, nPreRank + 1);
PreTotalSumW = Array{Float64}(nNullSimPts, nPreRank + 1);
tmpSumVec = Array{Float64}(nNullSimPts);
p1 = pointer(PrePartialSumW);
BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn), nNullSimPts), 1, p1, 1);
for i = 1 : nNullSimPts
  #pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
  #tmpSumVec[i] = BLAS.asum(1, pW, 1);
  tmpSumVec[i] = 0.0;
  PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
end
for j = 1 : nPreRank
  p1 = pointer(PrePartialSumW) + j * nNullSimPts * sizeof(Float64)
  BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn - j), nNullSimPts), 1, p1, 1);
  for i = 1 : nNullSimPts
    tmpSumVec[i] += WPreSim[j, i];
    PreTotalSumW[i, j+1] = tmpSumVec[i] + PrePartialSumW[i, j+1];
  end
end



if longitudinal == false && fine == false
     if test == "eRLRT"
        (b, vc0List, vc1List, pvalList) =
         vctest(ynew, [], XtNullBasis' * K , tests = test,
                Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                devices = device, PrePartialSumW = PrePartialSumW,
                PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                totalSumWConst = totalSumWConst, windowSize = tmpn - 1,
                partialSumW = partialSumW, totalSumW = totalSumW,
                lambda = lambda, W = W, nPreRank = nPreRank, nPtsChi2=nPtsChi2,
                offset = offset);

     else
        (b, vc0List, vc1List, pvalList) =
         vctest(y, X, K,tests = test, WPreSim = WPreSim,
                Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                devices = device, PrePartialSumW = PrePartialSumW,
                PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                totalSumWConst = totalSumWConst, windowSize = tmpn - 1,
                partialSumW = partialSumW, totalSumW = totalSumW,
                lambda = lambda, W = W, nPreRank = nPreRank, nPtsChi2=nPtsChi2,
                offset = offset);
     end
elseif longitudinal || fine
        (b, vc0List, vc1List, pvalList) =
         vctest(yWork, [], VWorkSqrt, tests = test, WPreSim = WPreSim,
                Vform = "half", nBlockAscent = nBlockAscent, nMMmax = nMMmax,
                nNullSimPts = nNullSimPts, pvalueComputings = pvalueComputing,
                nNullSimNewtonIter = nNullSimNewtonIter, tolX = tolX,
                devices = device, PrePartialSumW = PrePartialSumW,
                PreTotalSumW = PreTotalSumW, partialSumWConst = partialSumWConst,
                totalSumWConst = totalSumWConst, windowSize = tmpn - 1,
                partialSumW = partialSumW, totalSumW = totalSumW,
                lambda = lambda, W = W, nPreRank = nPreRank, nPtsChi2=nPtsChi2,
                offset = offset);
end

if out == true
#  outresults = vcat(vc0List, vc1List, pvalList)
  outrlst = hcat(yname, vcat(vc0List, vc1List, pvalList)')
  cname = hcat("phenotye", "vc0", "vc1", "pvalue")
  outrlst = [cname; outrlst]
  writedlm(outFile, outrlst)
end

return b, vc0List, vc1List, pvalList
end
