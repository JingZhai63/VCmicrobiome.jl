
#------------------------------------------------------------#
function MMicrobiomeVCTest(args...; 
                nObs::Int=86,
                device::String = "CPU",
                MemoryLimit::Int = 200000000, 
                covFile::String = "", 
                kernel::String = "multi",
                kernelFile::String = "",
                KernelLists::String = "",
                responseFile::String = "",
                yInit::Int = 5,
               
                obsN::Int=99,  # number of observations

                ZtZ::String = "intercept", 
                nMMmax::Int = 0,
                nBlockAscent::Int = 1000, # max block ascent iterations, default is 100
                nNullSimPts::Int = 10000, # simulation samples
                nNullSimNewtonIter::Int = 15, #Max Newton iteration
                test::String = "eLRT",
                outFile::String = "",
                tolX::Float64 = 1e-4,
                vcInit::Array{Float64, 1} = Float64[],
                Vform::String = "half",
                pvalueComputing::String = "chi2",pvalueComputings::String = "chi2",
                WPreSim::Array{Float64, 2} = [Float64[] Float64[]],
                PrePartialSumW::Array{Float64, 2} = [Float64[] Float64[]],
                PreTotalSumW::Array{Float64, 2} = [Float64[] Float64[]],
                partialSumWConst::Array{Float64, 1} = Float64[],
                totalSumWConst::Array{Float64, 1} = Float64[],
                windowSize::Int = 50,
                partialSumW::Array{Float64, 1} = Float64[],
                totalSumW::Array{Float64, 1} = Float64[],
                lambda::Array{Float64, 2} = [Float64[] Float64[]],
                W::Array{Float64, 2} = [Float64[] Float64[]],
                nPreRank::Int = 20,
                tmpmat0::Array{Float64, 2} = [Float64[] Float64[]],
                tmpmat1::Array{Float64, 2} = [Float64[] Float64[]],
                tmpmat2::Array{Float64, 2} = [Float64[] Float64[]],
                tmpmat3::Array{Float64, 2} = [Float64[] Float64[]],
                tmpmat4::Array{Float64, 2} = [Float64[] Float64[]],
                tmpmat5::Array{Float64, 2} = [Float64[] Float64[]],
                denomvec::Array{Float64, 1} = Float64[],
                d1f::Array{Float64, 1} = Float64[],
                d2f::Array{Float64, 1} = Float64[], offset::Int = 0,
                nPtsChi2::Int = 300,
                simnull::Vector{Float64} = Float64[])
 
  ##zj::-----------------------get kernel------------------------------#
  
  if kernel == "single"
    if isempty(kernelFile)
      error("mVctest:kernelwrongn\n",
               "# no kernel matrix")
    else
       K = readdlm(kernelFile,',');
       K=K[2:end,:];
       K=convert(Array{Float64, 2}, K);
       nObs=length(K[:,1]); # number of total observations N #58
    end
  elseif kernel == "multi"
       KK=KernlInput(KernelList=KernelLists,nObs=nObs);
       nObs =length(KK[1][:,1]);
  end

#zj::------------------ get covariate and response---------------------#
  if isempty(covFile)
    X = ones(nObs, 1); # N-by-1 matrix with all elements equal to 1
    X = convert(Array{Float64, 2}, X)
    if isempty(responseFile)
      error("MicrobiomVctest:noresponse\n", "# need to provide response file");
    else
      YY = readdlm(responseFile,',');
      YY = YY[2:end,yInit:end];  
      if size(YY, 1) != nObs
        error("mvctest:ywrongn\n",
              "# individuals in response file does not match dimension of V ");
      end
    end
  else
    X = readdlm(covFile,',');

    X = X[2:end,3:end];
    if
     size(X, 1) != nObs
      error("MicrobiomVctest:covwrongn\n",
            "# observationss in covariate file does not match kernel matrix");
    end
    if isempty(responseFile)
      error("MicrobiomVctest:noresponse\n",
            "# need to provide response file");
    else
      YY = readdlm(responseFile,',');
      YY = YY[2:end,yInit:end];
      
      if size(YY, 1) != nObs
        error("mvctest:ywrongn\n",
              "# individuals in response file does not match dimension of V ");
      end
    end

    X = convert(Array{Float64, 2}, X); # convert the Array type of X to Float64 with 2 dimentions
  end 
   simN=size(YY,2);

for iSim = 1:simN 
 
     y = YY[:,iSim];
     y = convert(Array{Float64, 1}, y); 
     nPerKeep = length(y);
     if isempty(covFile) # no covariates provided
       X = ones(nObs, 1); # N-by-1 matrix with all elements equal to 1
     else 
      X = readdlm(covFile,',');
      X = X[2:end,3:end];
     end
      X = convert(Array{Float64, 2}, X);
  nPerKeep=length(y);

    ## pre-processing for testing 

 if ZtZ == "none"
  # pre-compute basis of N(X') and projected y
    Xsvd = svdfact(X, thin = false);
    rankX = countnz(Xsvd[:S] .> nPerKeep * eps(Xsvd[:S][1]));
    XtNullBasis = Xsvd[:U][:, rankX + 1 : end]';
    ynew = zeros(nPerKeep - rankX);
    BLAS.gemv!('N', 1.0, XtNullBasis, y, 0.0, ynew);
    
    evalPhiAdj = Float64[];
    yShift = Float64[];
    XPhitNullBasis = [Float64[] Float64[]];
    rankQZ=0
    KPhiAdj = [Float64[] Float64[]];
    weightedW = [Float64[] Float64[]];
    QZ=[Float64[] Float64[]];
    nPreRank=nPerKeep; #?????
    PrePartialSumW = [Float64[] Float64[]];
    PreTotalSumW = [Float64[] Float64[]];
    tmpn = nObs - rankX; ##zj::

    partialSumWConst = Array(Float64, nNullSimPts);
    totalSumWConst = Array(Float64, nNullSimPts);
    offset = 0;
 elseif ZtZ =="intercept"
       # >=2 non-trivial variance components
     ynew = Float64[];
     Covariates= readdlm(covFile,',');
     CovariatesID = Covariates[2:end,1];
     Nsort=sort(unique(CovariatesID));

     counts = map(x->count(y->x==y,CovariatesID),Nsort);

     Z=Array(Float64,length(CovariatesID),length(counts));
     fill!(Z,0.0);
     Count=counts[1];
     Z[1:counts[1],1]=ones(counts[1]);
     for i = 2:length(counts)
      Z[(Count+1):(Count+counts[i]),i]=ones(counts[i]);
      Count=Count+counts[i];
     end

    # obtain a basis of X and (possibly) extra variance components
    # Xsvd[:U] -- QX
    Xsvd = svdfact(X, thin = false);
    rankX = countnz(Xsvd[:S] .> nPerKeep * eps(Xsvd[:S][1])); # r0=rank(Q0)
    XtNullBasis = Xsvd[:U][:, rankX + 1 : end]; 
    QX = Xsvd[:U][:, 1 : rankX];  #orthonormal basis of Q0 of C(X)

    
 QZsvd = svdfact(BLAS.gemm('T', 'N', XtNullBasis, Z)); # XtNullBasis N(X')
    ZKeepIdx = (1 - QZsvd[:S]) .< 1e-6;
    rankQZ = countnz(ZKeepIdx);
    QZ = Array(Float64, nPerKeep, rankQZ); 
    BLAS.gemm!('N', 'N', 1.0, Z, QZsvd[:V][:, 1 : rankQZ], 0.0, QZ); #orthonormal basis Q1 of C(Phi)-C(X)
    XPhitNullBasis = null([QX Z]');  

    tmpMat = Array(Float64, rankQZ, size(Z,2)); 
    BLAS.gemm!('T', 'N', 1.0, QZ, Z, 0.0, tmpMat); 
    PhiAdjsvd = svdfact(tmpMat, thin = false); #PhiAdjsvd ---WLambdaW'
   
    #W = PhiAdjsvd[:U];
    evalPhiAdj = PhiAdjsvd[:S] .^ 2;
    # enforce first entry of each eigenvector to be >=0
    idxW = vec(PhiAdjsvd[:U][1, :]) .< 0;
    PhiAdjsvd[:U][:, idxW] = - PhiAdjsvd[:U][:, idxW];
    #scale!(PhiAdjsvd[:U][:, idxW], -1.0);
    #KPhiAdj = PhiAdjsvd[:U] .* sqrt(evalPhiAdj' / minimum(evalPhiAdj) - 1);
    KPhiAdj = PhiAdjsvd[:U][:, :]; # KPhiAdj-- W  95 95
    scale!(KPhiAdj, sqrt(evalPhiAdj / minimum(evalPhiAdj) - 1));  ##KPhiAdj---K
    InvSqrtevalPhiAdj = similar(evalPhiAdj);
    for i = 1 : length(evalPhiAdj)
      InvSqrtevalPhiAdj[i] = 1 / sqrt(evalPhiAdj[i]);
    end
    weightedW = PhiAdjsvd[:U][:, :];
    scale!(weightedW, InvSqrtevalPhiAdj);  ##weightedW---W'*Lambda^1/2

    # precompute shift in Y
    yShift = QZ' * y;

    # prepare for simulating Chi Squares and sums of Chi Squares
    tmpn = length(evalPhiAdj);
  end

  if isempty(outFile)
    outFile = string(responseFile,test,iSim,".out");
  else
    outFile = string(outFile,".out")
  end
  fid = open(outFile, "w");
  println(fid, "vc0,vc1,Pvalue");
  close(fid);

#### test group by group

#    outresults = loopKernel()
   if ZtZ == "intercept"
     QRes = Array(Float64, nPerKeep, rankQZ); ## Q2

     tmpvec = similar(yShift);
     tmpvecQRes = Array(Float64, rankQZ); ##
     yWork = similar(evalPhiAdj);
     partialSumWConst = Array(Float64, nNullSimPts);
     totalSumWConst = Array(Float64, nNullSimPts);
     subXPhitSV = Array(Float64, size(XPhitNullBasis, 2), rankQZ);
     pSubXPhitSV = pointer(subXPhitSV);
     offset = 0;

     tmpMat = Array(Float64, nPerKeep, size(XPhitNullBasis, 2)); #windowSize replaced by nPerKeep
     
     ##change windowSize to nPerKeep
     VWorkSqrt = Array(Float64, length(evalPhiAdj), nPerKeep); 
     VWorkSqrt2 = Array(Float64, length(evalPhiAdj), nPerKeep);
   end
 
  if kernel == "multi"
    nGrps = length(KK);
  elseif kernel == "single"
    nGrps = 1;
  end
  
  results = Array(Any, nGrps, 3);
  idx = 0;

    vc0List = zeros(nGrps);
    vc1List = zeros(nGrps);
    pvalList = zeros(nGrps);

for g = 1:nGrps
     y = YY[:,iSim];
     y = convert(Array{Float64, 1}, y); 
     nPerKeep = length(y);
  
   if kernel == "multi"
     K=KK[g];
   end

  if ZtZ == "none"

        if test == "eRLRT"
              #----------------testvc1 is for eRLRT (baseline)---------------#    
    function testvc1(y=ynew,V=XtNullBasis*K, X=[],ZtZ=ZtZ,tests=test, devices=device, vcInit = Float64[],
                       bInit = Float64[], WPreSim= [Float64[] Float64[]],
                       pvalueComputings = pvalueComputing,nMMmax=nMMmax)
       
         n = length(y);
         if size(V, 1) != n
           error("vctest:wrongdimV\n", "dimension of V does not match that of X");
         end
       
         # set default maximum MM iteration
         if nMMmax == 0 && tests == "eLRT"
           nMMmax = 10;
         elseif nMMmax == 0 && tests == "eRLRT"
           nMMmax = 1000;
         end
       
         # SVD of X
         if isempty(X)
           rankX = 0;
           X = reshape(X, n, 0);
           X = convert(Array{Float64, 2}, X);
           # LRT is same as RLRT if X is null
           if tests == "eLRT"
             tests = "eRLRT";
           end
         else
           if tests == "eRLRT"
             (UX, svalX) = svd(X, thin = false);
           else
             (UX, svalX) = svd(X);
           end
           rankX = countnz(svalX .> n * eps(svalX[1]));
         end
       
         # eigendecomposition of V
         if Vform == "whole"
           (evalV, UV) = eig(V);
           rankV = countnz(evalV .> n * eps(sqrt(maximum(evalV))));
           sortIdx = sortperm(evalV, rev = true);
           evalV = evalV[sortIdx[1:rankV]];
           UV = UV[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = V;
           end
         elseif Vform == "half"
           (UVfull, tmpevalVfull) = svd(V, thin=false);
           evalVfull = zeros(n);
           pevalVfull = pointer(evalVfull);
           ptmpevalV = pointer(tmpevalVfull);
           BLAS.blascopy!(length(tmpevalVfull), ptmpevalV, 1, pevalVfull, 1);
           rankV = countnz(evalVfull .> n * eps(maximum(evalVfull)));
           evalVfull = evalVfull .^ 2;
           evalV = evalVfull[1:rankV];
           UV = UVfull[:, 1:rankV];
           #if tests == "eLRT" || tests == "eRLRT"
           #  wholeV = *(V, V');
           #end
         elseif Vform == "eigen"
           UV = V.U;
           evalV = V.eval;
           rankV = countnz(V.eval .> n * eps(sqrt(maximum(V.eval))));
           sortIdx = sortperm(V.eval, rev = true);
           V.eval = V.eval[sortIdx[1:rankV]];
           V.U = V.U[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = *(V.U .* reshape(V.eval, length(V.eval), 1), V.U');
           end
         end
       
         # obtain eigenvalues of (I-PX)V(I-PX)
         if !isempty(X) || tests == "eScore"
           #sqrtV = UV .* sqrt(evalV)';
           sqrtV = similar(UV);
           psqrtV = pointer(sqrtV);
           pUV = pointer(UV);
           BLAS.blascopy!(n*rankV, pUV, 1, psqrtV, 1);
           scale!(sqrtV, sqrt(evalV));
           # scale!(sqrtV, sqrt(abs(evalV)));
         end
         if isempty(X)
           evalAdjV = evalV;
         else
           subUX = Array(Float64, n, rankX);
           psubUX = pointer(subUX);
           pUX = pointer(UX);
           BLAS.blascopy!(n*rankX, pUX, 1, psubUX, 1);
           mat1 = BLAS.gemm('T', 'N', 1.0, subUX, sqrtV);
           mat2 = BLAS.gemm('N', 'N', 1.0, subUX, mat1);
           (UAdjV, evalAdjV) = svd(sqrtV - mat2, thin = false);
           if isempty(evalAdjV)
             evalAdjV = Float64[];
           else
             evalAdjV = evalAdjV[evalAdjV .> n * eps(maximum(evalAdjV))] .^ 2;
           end
         end
         rankAdjV = length(evalAdjV);
       
         nSimPts = nNullSimPts;
         if size(WPreSim, 1) < rankAdjV
             newSim = randn(rankAdjV - size(WPreSim, 1), nSimPts) .^ 2;
             if isempty(WPreSim)
               WPreSim = newSim;    ##WPreSim=54,10000
             else
               WPreSim = [WPreSim, newSim];
             end
             #windowSize = rankAdjV;
         end
       
         if tests == "eRLRT" || tests == "eScore"
          windowSize=rankAdjV;
         elseif tests == "eLRT"
          windowSize=rankV;
         end
       
         if ZtZ =="intercept"
         nPreRank = min(length(evalPhiAdj), windowSize);
         else
          nPreRank=n-rankX;  
         end
        
         PrePartialSumW = Array(Float64, nNullSimPts, nPreRank+1);
         PreTotalSumW = Array(Float64, nNullSimPts, nPreRank+1);
         tmpSumVec = Array(Float64, nNullSimPts);
         p1 = pointer(PrePartialSumW);
         BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn), nNullSimPts), 1, p1, 1);
         for i = 1 : nNullSimPts
           #pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
           #tmpSumVec[i] = BLAS.asum(1, pW, 1);
           tmpSumVec[i] = 0.0;
           PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
         end
         for j = 1 : nPreRank
           p1 = pointer(PrePartialSumW) + j * nNullSimPts * sizeof(Float64);
           BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn +1 - j), nNullSimPts), 1, p1, 1); ##ZJ:: tmpn-j???
           for i = 1 : nNullSimPts
             tmpSumVec[i] += WPreSim[j, i];
             PreTotalSumW[i, j+1] = tmpSumVec[i] + PrePartialSumW[i, j+1];
           end
         end
       
         if pvalueComputing == "chi2"
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nPtsChi2);
           totalSumW = Array(Float64, nPtsChi2);
           lambda = Array(Float64, 1, nPtsChi2);
           W = Array(Float64, windowSize, nPtsChi2);
           tmpmat0 = Array(Float64, windowSize, nPtsChi2);
           tmpmat1 = Array(Float64, windowSize, nPtsChi2);
           tmpmat2 = Array(Float64, windowSize, nPtsChi2);
           tmpmat3 = Array(Float64, windowSize, nPtsChi2);
           tmpmat4 = Array(Float64, windowSize, nPtsChi2);
           tmpmat5 = Array(Float64, windowSize, nPtsChi2);
           denomvec = Array(Float64, nPtsChi2);
           d1f = Array(Float64, nPtsChi2);
           d2f = Array(Float64, nPtsChi2);
           if test == "eScore"
             simnull = Array(Float64, nPtsChi2);
           else
             simnull = Float64[];
           end
         else
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nNullSimPts);
           totalSumW = Array(Float64, nNullSimPts);
           lambda = Array(Float64, 1, nNullSimPts);
           W = Array(Float64, windowSize, nNullSimPts);
           tmpmat0 = Array(Float64, windowSize, nNullSimPts);
           tmpmat1 = Array(Float64, windowSize, nNullSimPts);
           tmpmat2 = Array(Float64, windowSize, nNullSimPts);
           tmpmat3 = Array(Float64, windowSize, nNullSimPts);
           tmpmat4 = Array(Float64, windowSize, nNullSimPts);
           tmpmat5 = Array(Float64, windowSize, nNullSimPts);
           denomvec = Array(Float64, nNullSimPts);
           d1f = Array(Float64, nNullSimPts);
           d2f = Array(Float64, nNullSimPts);
           if test == "eScore"
             simnull = Array(Float64, nNullSimPts);
           else
             simnull = Float64[];
           end
         end
       #end
       
         # fit the variance component model
         if tests == "eLRT"
       
           # estimates under null model
           bNull = X \ y;
           rNull = y - X * bNull;
           vc0Null = norm(rNull) ^ 2 / n;
           Xrot = UVfull' * X;
           yrot = UVfull' * y;
           loglConst = - 0.5 * n * log(2.0 * pi);
       
           # set starting point
           if isempty(bInit)
             b = copy(bNull);
             r = copy(rNull);
           else
             b = copy(bInit);
             r = y - X * b;
           end
           if isempty(vcInit)
             vc0 = norm(r) ^ 2 / n;
             vc1 = 1;
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
           else
             vc0 = vcInit[1];
             # avoid sticky boundary
             if vc0 == 0
               vc0 = 1e-4;
             end
             vc1 = vcInit[2];
             if vc1 == 0
               vc1 = 1e-4;
             end
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
           end
       
           # update residuals according supplied var. components
           r = y - BLAS.gemv('N', X, b);
           rnew = BLAS.gemv('T', UVfull, r);
           #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
           logl = loglConst + sum(log, wt) - 0.5 * sumabs2(rnew .* wt);
       
           nIters = 0;
           #bOld = similar(b);
           #pbOld = pointer(bOld);
           #pb = pointer(b);
           #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
           #denvec = similar(rnew);
           denvec = Float64[];
           for iBlockAscent = 1:nBlockAscent
       
             nIters = iBlockAscent;
       
             # update variance components
             for iMM = 1:nMMmax
               denvec = 1.0 ./ (vc1 * evalVfull + vc0);
               numvec = rnew .* denvec;
               #vc0 = sqrt( (vc0 ^ 2 * sumabs2(numvec) + deltaRes) /
               #            (sumabs(denvec) + (n - rankV) / vc0) );
               vc0 = vc0 * sqrt(sumabs2(numvec) / sumabs(denvec));
               vc1 = vc1 * sqrt( dot(numvec, numvec .* evalVfull) /
                                  sumabs(evalVfull .* denvec) );
               #wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             end
       
             # update fixed effects and residuals
             #pb = pointer(b);
             #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
             wt = sqrt(denvec);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
             r = y - BLAS.gemv('N', X, b);
             rnew = BLAS.gemv('T', UVfull, r);
             #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
       
             # stopping criterion
             loglOld = logl;
             logl = loglConst + sum(log, wt) -
               0.5 * BLAS.dot(length(evalVfull), rnew .^ 2, 1, denvec, 1);
             if abs(logl - loglOld) < tolX * (abs(logl) + 1.0)
               break
             end
       
           end
       
           # log-likelihood at alternative
           logLikeAlt = logl;
           # log-likelihood at null
           logLikeNull = loglConst - 0.5 * n * log(vc0Null) -
             0.5 * sum(rNull .^ 2) / vc0Null;
           if logLikeNull >= logLikeAlt
             vc0 = vc0Null;
             vc1 = 0;
             logLikeAlt = logLikeNull;
             b = bNull;
             r = rNull;
           end
       
           # LRT test statistic
           statLRT = 2 * (logLikeAlt - logLikeNull);
       
           # obtain p-value for testing vc1=0
           vc1_pvalue = vctestnullsim(statLRT, evalV, evalAdjV, n, rankX,
                                      WPreSim, device = devices,
                                      nSimPts = nNullSimPts,
                                      nNewtonIter = nNullSimNewtonIter,
                                      test = "eLRT",
                                      pvalueComputing = pvalueComputings,
                                      PrePartialSumW = PrePartialSumW,
                                      PreTotalSumW = PreTotalSumW,
                                      partialSumWConst = partialSumWConst,
                                      totalSumWConst = totalSumWConst,
                                      windowSize = windowSize,
                                      partialSumW = partialSumW,
                                      totalSumW = totalSumW,
                                      lambda = lambda, W = W,
                                      nPreRank = nPreRank,
                                      tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                      tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                      tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                      denomvec = denomvec,
                                      d1f = d1f, d2f = d2f, offset = offset,
                                      nPtsChi2 = nPtsChi2, simnull = simnull);
       
           # return values
           return b, vc0, vc1, vc1_pvalue;
       
         elseif tests == "eRLRT"
       
           if isempty(X)
             ytilde = y;
             rankBVB = rankV;
             evalBVB = evalV;
             UBVB = UV;
           else
             # obtain a basis of N(X')
             B = UX[:, rankX+1:end];
             ytilde = B' * y;
       
             # eigen-decomposition of B'VB and transform data
             (UBVB, evalBVB) = svd(B' * sqrtV);
             rankBVB = countnz(evalBVB .> n * eps(maximum(evalBVB)));
             evalBVB = evalBVB[1:rankBVB] .^ 2;
             UBVB = UBVB[:, 1:rankBVB];
           end
           resnew = UBVB' * ytilde;
           normYtilde2 = norm(ytilde) ^ 2;
           deltaRes = normYtilde2 - norm(resnew) ^ 2;
       
           # set initial point
           # TODO: is there better initial point?
           vc0Null = normYtilde2 / length(ytilde);
           if isempty(vcInit)
             vc0 = vc0Null;
             vc1 = 1.0;
           else
             vc0 = vcInit[1];
             vc1 = vcInit[2];
           end
       
      # MM loop for estimating variance components
            nIters = 0;
            #tmpvec = Array(Float64, rankBVB);
            #denvec = Array(Float64, rankBVB);
            #numvec = Array(Float64, rankBVB);
            tmpvec = Float64[];
            for iMM = 1:nMMmax
              nIters = iMM;
              #tmpvec = vc0 + vc1 * evalBVB;
              #tmpvec = evalBVB;
              #numvecSum = 0.0;
              #denvecSum = 0.0;
              #numvecProdSum = 0.0;
              #denvecProdSum = 0.0;
              tmpvec = BLAS.scal(rankBVB, vc1, evalBVB, 1);
              #for i = 1 : rankBVB
              #  tmpvec[i] += vc0;
              #  denvec[i] = 1 / tmpvec[i];
              #  numvec[i] = (resnew[i] * denvec[i]) ^ 2;
              #  numvecSum += numvec[i];
              #  denvecSum += denvec[i];
              #  numvecProdSum += evalBVB[i] * numvec[i];
              #  denvecProdSum += evalBVB[i] * denvec[i];
              #end
              tmpvec += vc0;
              denvec = 1 ./ tmpvec;
              numvec = (resnew .* denvec) .^ 2;
              vcOld = [vc0 vc1];
              vc0 = sqrt( (vc0 ^ 2 * sum(numvec) + deltaRes) /
                           (sum(denvec) + (n - rankX - rankBVB) / vc0) );
              #vc0 = sqrt( (vc0 ^ 2 * numvecSum + deltaRes) /
              #             (denvecSum + (n - rankX - rankBVB) / vc0) );
              vc1 = vc1 * sqrt( sum(evalBVB .* numvec) / sum(evalBVB .* denvec) );
              #vc1 = vc1 * sqrt( numvecProdSum / denvecProdSum );
              # stopping criterion
              if norm([vc0 vc1] - vcOld) <= tolX * (norm(vcOld) + 1)
                break;
              end
            end
        
            # restrictive log-likelihood at alternative
        
            loglConst = - 0.5 * (n - rankX) * log(2 * pi);
            logLikeAlt =  loglConst - 0.5 * sum(log(tmpvec)) -
              0.5 * (n - rankX - rankBVB) * log(vc0) - 0.5 * normYtilde2 / vc0 +
              0.5 * sum(resnew .^ 2 .* (1 / vc0 - 1 ./ (tmpvec)));
            # restrictive log-likelihood at null
            logLikeNull = - 0.5 * (n - rankX) * (log(2 * pi) + log(vc0Null)) -
              0.5 / vc0Null * normYtilde2;
            if logLikeNull >= logLikeAlt
              vc0 = vc0Null;
              vc1 = 0;
              logLikeAlt = logLikeNull;
            end
        
            # RLRT test statitic
            statRLRT = 2 * (logLikeAlt - logLikeNull);
        
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statRLRT, evalV, evalAdjV, n, rankX,
                                       WPreSim, device = devices,
                                       nSimPts = nNullSimPts,
                                       nNewtonIter = nNullSimNewtonIter,
                                       test = "eRLRT",
                                       pvalueComputing = pvalueComputings,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
        
            # estimate fixed effects
            if isempty(X)
              b = zeros(0);
            else
              Xrot = UV' * X;
              yrot = UV' * y;
              wt = 1.0 ./ sqrt(vc1 * evalV + vc0);
              Xnew = scale(wt, Xrot);
              ynew = wt .* yrot;
              b = Xnew \ ynew;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          elseif tests == "eScore"
        
            # fit the null model
            b = X \ y;
            r = y - X * b;
            vc0 = norm(r) ^ 2 / n;
            vc1 = 0;
        
            # score test statistic
            #statScore = norm(r' * sqrtV) ^ 2 / norm(r) ^ 2;
            statScore = norm(BLAS.gemv('T', sqrtV, r)) ^ 2 / norm(r) ^ 2;
            #statScore = max(statScore, sum(evalV) / n);
        
            #=
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statScore, evalV, evalAdjV, n, rankX,
                                       WPreSim, test = "eScore",
                                       nSimPts = nNullSimPts,
                                       pvalueComputing = pvalueComputings,
                                       nNewtonIter = nNullSimNewtonIter,
                                       device = devices,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
            =#
        
            if statScore <= sum(evalV) / n
              vc1_pvalue = 1.0;
            else
              w = [evalAdjV - statScore; - statScore];
              dfvec = [ones(Int, rankAdjV); n - rankX - rankAdjV];
              fun(x) = imag(prod((1 - 2 .* w * x * im) .^ (- dfvec / 2)) / x);
              (vc1_pvalue, err) = quadgk(fun, 0, Inf);
              vc1_pvalue = 0.5 + vc1_pvalue / pi;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          end
        end
 #------- 
          (b, vc0List[g], vc1List[g], pvalList[g]) = testvc1();
        else
          #----------------testvc2 is for eLRT and eScore(baseline)---------------#    
    

    function testvc2(y=y,V=K, X=X,tests=test, devices=device, vcInit = Float64[],
                       bInit = Float64[], WPreSim= [Float64[] Float64[]],
                       pvalueComputings = pvalueComputing,nMMmax=nMMmax)
       
         n = length(y);
         if size(V, 1) != n
           error("vctest:wrongdimV\n", "dimension of V does not match that of X");
         end
       
         # set default maximum MM iteration
         if nMMmax == 0 && tests == "eLRT"
           nMMmax = 10;
         elseif nMMmax == 0 && tests == "eRLRT"
           nMMmax = 1000;
         end
       
         # SVD of X
         if isempty(X)
           rankX = 0;
           X = reshape(X, n, 0);
           X = convert(Array{Float64, 2}, X);
           # LRT is same as RLRT if X is null
           if tests == "eLRT"
             tests = "eRLRT";
           end
         else
           if tests == "eRLRT"
             (UX, svalX) = svd(X, thin = false);
           else
             (UX, svalX) = svd(X);
           end
           rankX = countnz(svalX .> n * eps(svalX[1]));
         end
       
         # eigendecomposition of V
         if Vform == "whole"
           (evalV, UV) = eig(V);
           rankV = countnz(evalV .> n * eps(sqrt(maximum(evalV))));
           sortIdx = sortperm(evalV, rev = true);
           evalV = evalV[sortIdx[1:rankV]];
           UV = UV[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = V;
           end
         elseif Vform == "half"
           (UVfull, tmpevalVfull) = svd(V, thin=false);
           evalVfull = zeros(n);
           pevalVfull = pointer(evalVfull);
           ptmpevalV = pointer(tmpevalVfull);
           BLAS.blascopy!(length(tmpevalVfull), ptmpevalV, 1, pevalVfull, 1);
           rankV = countnz(evalVfull .> n * eps(maximum(evalVfull)));
           evalVfull = evalVfull .^ 2;
           evalV = evalVfull[1:rankV];
           UV = UVfull[:, 1:rankV];
           #if tests == "eLRT" || tests == "eRLRT"
           #  wholeV = *(V, V');
           #end
         elseif Vform == "eigen"
           UV = V.U;
           evalV = V.eval;
           rankV = countnz(V.eval .> n * eps(sqrt(maximum(V.eval))));
           sortIdx = sortperm(V.eval, rev = true);
           V.eval = V.eval[sortIdx[1:rankV]];
           V.U = V.U[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = *(V.U .* reshape(V.eval, length(V.eval), 1), V.U');
           end
         end
       
         # obtain eigenvalues of (I-PX)V(I-PX)
         if !isempty(X) || tests == "eScore"
           #sqrtV = UV .* sqrt(evalV)';
           sqrtV = similar(UV);
           psqrtV = pointer(sqrtV);
           pUV = pointer(UV);
           BLAS.blascopy!(n*rankV, pUV, 1, psqrtV, 1);
           scale!(sqrtV, sqrt(evalV));
           # scale!(sqrtV, sqrt(abs(evalV)));
         end
         if isempty(X)
           evalAdjV = evalV;
         else
           subUX = Array(Float64, n, rankX);
           psubUX = pointer(subUX);
           pUX = pointer(UX);
           BLAS.blascopy!(n*rankX, pUX, 1, psubUX, 1);
           mat1 = BLAS.gemm('T', 'N', 1.0, subUX, sqrtV);
           mat2 = BLAS.gemm('N', 'N', 1.0, subUX, mat1);
           (UAdjV, evalAdjV) = svd(sqrtV - mat2, thin = false);
           if isempty(evalAdjV)
             evalAdjV = Float64[];
           else
             evalAdjV = evalAdjV[evalAdjV .> n * eps(maximum(evalAdjV))] .^ 2;
           end
         end
         rankAdjV = length(evalAdjV);
       
         srand(1)
         nSimPts = nNullSimPts;
         if size(WPreSim, 1) < rankAdjV
             newSim = randn(rankAdjV - size(WPreSim, 1), nSimPts) .^ 2;
             if isempty(WPreSim)
               WPreSim = newSim;    ##WPreSim=54,10000
             else
               WPreSim = [WPreSim, newSim];
             end
             #windowSize = rankAdjV;
         end
       
         if test == "eRLRT" || test == "eScore"
          windowSize=rankAdjV;
         elseif test == "eLRT"
          windowSize=rankV;
         end
       
         if ZtZ =="intercept"
         nPreRank = min(length(evalPhiAdj), windowSize);
         else
          nPreRank=n-rankX;  #zj::
         end
         # nPreRank= rankAdjV
        
         PrePartialSumW = Array(Float64, nNullSimPts, nPreRank+1);
         PreTotalSumW = Array(Float64, nNullSimPts, nPreRank+1);
         tmpSumVec = Array(Float64, nNullSimPts);
         p1 = pointer(PrePartialSumW);
         BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn), nNullSimPts), 1, p1, 1);
         for i = 1 : nNullSimPts
           #pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
           #tmpSumVec[i] = BLAS.asum(1, pW, 1);
           tmpSumVec[i] = 0.0;
           PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
         end
         for j = 1 : nPreRank
           p1 = pointer(PrePartialSumW) + j * nNullSimPts * sizeof(Float64);
           BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn +1 - j), nNullSimPts), 1, p1, 1); ##ZJ:: tmpn-j???
           for i = 1 : nNullSimPts
             tmpSumVec[i] += WPreSim[j, i];
             PreTotalSumW[i, j+1] = tmpSumVec[i] + PrePartialSumW[i, j+1];
           end
         end
       
         if pvalueComputing == "chi2"
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nPtsChi2);
           totalSumW = Array(Float64, nPtsChi2);
           lambda = Array(Float64, 1, nPtsChi2);
           W = Array(Float64, windowSize, nPtsChi2);
           tmpmat0 = Array(Float64, windowSize, nPtsChi2);
           tmpmat1 = Array(Float64, windowSize, nPtsChi2);
           tmpmat2 = Array(Float64, windowSize, nPtsChi2);
           tmpmat3 = Array(Float64, windowSize, nPtsChi2);
           tmpmat4 = Array(Float64, windowSize, nPtsChi2);
           tmpmat5 = Array(Float64, windowSize, nPtsChi2);
           denomvec = Array(Float64, nPtsChi2);
           d1f = Array(Float64, nPtsChi2);
           d2f = Array(Float64, nPtsChi2);
           if test == "eScore"
             simnull = Array(Float64, nPtsChi2);
           else
             simnull = Float64[];
           end
         else
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nNullSimPts);
           totalSumW = Array(Float64, nNullSimPts);
           lambda = Array(Float64, 1, nNullSimPts);
           W = Array(Float64, windowSize, nNullSimPts);
           tmpmat0 = Array(Float64, windowSize, nNullSimPts);
           tmpmat1 = Array(Float64, windowSize, nNullSimPts);
           tmpmat2 = Array(Float64, windowSize, nNullSimPts);
           tmpmat3 = Array(Float64, windowSize, nNullSimPts);
           tmpmat4 = Array(Float64, windowSize, nNullSimPts);
           tmpmat5 = Array(Float64, windowSize, nNullSimPts);
           denomvec = Array(Float64, nNullSimPts);
           d1f = Array(Float64, nNullSimPts);
           d2f = Array(Float64, nNullSimPts);
           if test == "eScore"
             simnull = Array(Float64, nNullSimPts);
           else
             simnull = Float64[];
           end
         end
       #end
       
         # fit the variance component model
         if tests == "eLRT"
       
           # estimates under null model
           bNull = X \ y;
           rNull = y - X * bNull;
           vc0Null = norm(rNull) ^ 2 / n;
           Xrot = UVfull' * X;
           yrot = UVfull' * y;
           loglConst = - 0.5 * n * log(2.0 * pi);
       
           # set starting point
           if isempty(bInit)
             b = copy(bNull);
             r = copy(rNull);
           else
             b = copy(bInit);
             r = y - X * b;
           end
           if isempty(vcInit)
             vc0 = norm(r) ^ 2 / n;
             vc1 = 1;
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
           else
             vc0 = vcInit[1];
             # avoid sticky boundary
             if vc0 == 0
               vc0 = 1e-4;
             end
             vc1 = vcInit[2];
             if vc1 == 0
               vc1 = 1e-4;
             end
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
           end
       
           # update residuals according supplied var. components
           r = y - BLAS.gemv('N', X, b);
           rnew = BLAS.gemv('T', UVfull, r);
           #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
           logl = loglConst + sum(log, wt) - 0.5 * sumabs2(rnew .* wt);
       
           nIters = 0;
           #bOld = similar(b);
           #pbOld = pointer(bOld);
           #pb = pointer(b);
           #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
           #denvec = similar(rnew);
           denvec = Float64[];
           for iBlockAscent = 1:nBlockAscent
       
             nIters = iBlockAscent;
       
             # update variance components
             for iMM = 1:nMMmax
               denvec = 1.0 ./ (vc1 * evalVfull + vc0);
               numvec = rnew .* denvec;
               #vc0 = sqrt( (vc0 ^ 2 * sumabs2(numvec) + deltaRes) /
               #            (sumabs(denvec) + (n - rankV) / vc0) );
               vc0 = vc0 * sqrt(sumabs2(numvec) / sumabs(denvec));
               vc1 = vc1 * sqrt( dot(numvec, numvec .* evalVfull) /
                                  sumabs(evalVfull .* denvec) );
               #wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             end
       
             # update fixed effects and residuals
             #pb = pointer(b);
             #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
             wt = sqrt(denvec);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
             r = y - BLAS.gemv('N', X, b);
             rnew = BLAS.gemv('T', UVfull, r);
             #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
       
             # stopping criterion
             loglOld = logl;
             logl = loglConst + sum(log, wt) -
               0.5 * BLAS.dot(length(evalVfull), rnew .^ 2, 1, denvec, 1);
             if abs(logl - loglOld) < tolX * (abs(logl) + 1.0)
               break
             end
       
           end
       
           # log-likelihood at alternative
           logLikeAlt = logl;
           # log-likelihood at null
           logLikeNull = loglConst - 0.5 * n * log(vc0Null) -
             0.5 * sum(rNull .^ 2) / vc0Null;
           if logLikeNull >= logLikeAlt
             vc0 = vc0Null;
             vc1 = 0;
             logLikeAlt = logLikeNull;
             b = bNull;
             r = rNull;
           end
       
           # LRT test statistic
           statLRT = 2 * (logLikeAlt - logLikeNull);
       
           # obtain p-value for testing vc1=0
           vc1_pvalue = vctestnullsim(statLRT, evalV, evalAdjV, n, rankX,
                                      WPreSim, device = devices,
                                      nSimPts = nNullSimPts,
                                      nNewtonIter = nNullSimNewtonIter,
                                      test = "eLRT",
                                      pvalueComputing = pvalueComputings,
                                      PrePartialSumW = PrePartialSumW,
                                      PreTotalSumW = PreTotalSumW,
                                      partialSumWConst = partialSumWConst,
                                      totalSumWConst = totalSumWConst,
                                      windowSize = windowSize,
                                      partialSumW = partialSumW,
                                      totalSumW = totalSumW,
                                      lambda = lambda, W = W,
                                      nPreRank = nPreRank,
                                      tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                      tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                      tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                      denomvec = denomvec,
                                      d1f = d1f, d2f = d2f, offset = offset,
                                      nPtsChi2 = nPtsChi2, simnull = simnull);
       
           # return values
           return b, vc0, vc1, vc1_pvalue;
       
         elseif tests == "eRLRT"
       
           if isempty(X)
             ytilde = y;
             rankBVB = rankV;
             evalBVB = evalV;
             UBVB = UV;
           else
             # obtain a basis of N(X')
             B = UX[:, rankX+1:end];
             ytilde = B' * y;
       
             # eigen-decomposition of B'VB and transform data
             (UBVB, evalBVB) = svd(B' * sqrtV);
             rankBVB = countnz(evalBVB .> n * eps(maximum(evalBVB)));
             evalBVB = evalBVB[1:rankBVB] .^ 2;
             UBVB = UBVB[:, 1:rankBVB];
           end
           resnew = UBVB' * ytilde;
           normYtilde2 = norm(ytilde) ^ 2;
           deltaRes = normYtilde2 - norm(resnew) ^ 2;
       
           # set initial point
           # TODO: is there better initial point?
           vc0Null = normYtilde2 / length(ytilde);
           if isempty(vcInit)
             vc0 = vc0Null;
             vc1 = 1.0;
           else
             vc0 = vcInit[1];
             vc1 = vcInit[2];
           end
       
      # MM loop for estimating variance components
            nIters = 0;
            #tmpvec = Array(Float64, rankBVB);
            #denvec = Array(Float64, rankBVB);
            #numvec = Array(Float64, rankBVB);
            tmpvec = Float64[];
            for iMM = 1:nMMmax
              nIters = iMM;
              #tmpvec = vc0 + vc1 * evalBVB;
              #tmpvec = evalBVB;
              #numvecSum = 0.0;
              #denvecSum = 0.0;
              #numvecProdSum = 0.0;
              #denvecProdSum = 0.0;
              tmpvec = BLAS.scal(rankBVB, vc1, evalBVB, 1);
              #for i = 1 : rankBVB
              #  tmpvec[i] += vc0;
              #  denvec[i] = 1 / tmpvec[i];
              #  numvec[i] = (resnew[i] * denvec[i]) ^ 2;
              #  numvecSum += numvec[i];
              #  denvecSum += denvec[i];
              #  numvecProdSum += evalBVB[i] * numvec[i];
              #  denvecProdSum += evalBVB[i] * denvec[i];
              #end
              tmpvec += vc0;
              denvec = 1 ./ tmpvec;
              numvec = (resnew .* denvec) .^ 2;
              vcOld = [vc0 vc1];
              vc0 = sqrt( (vc0 ^ 2 * sum(numvec) + deltaRes) /
                           (sum(denvec) + (n - rankX - rankBVB) / vc0) );
              #vc0 = sqrt( (vc0 ^ 2 * numvecSum + deltaRes) /
              #             (denvecSum + (n - rankX - rankBVB) / vc0) );
              vc1 = vc1 * sqrt( sum(evalBVB .* numvec) / sum(evalBVB .* denvec) );
              #vc1 = vc1 * sqrt( numvecProdSum / denvecProdSum );
              # stopping criterion
              if norm([vc0 vc1] - vcOld) <= tolX * (norm(vcOld) + 1)
                break;
              end
            end
        
            # restrictive log-likelihood at alternative
        
            loglConst = - 0.5 * (n - rankX) * log(2 * pi);
            logLikeAlt =  loglConst - 0.5 * sum(log(tmpvec)) -
              0.5 * (n - rankX - rankBVB) * log(vc0) - 0.5 * normYtilde2 / vc0 +
              0.5 * sum(resnew .^ 2 .* (1 / vc0 - 1 ./ (tmpvec)));
            # restrictive log-likelihood at null
            logLikeNull = - 0.5 * (n - rankX) * (log(2 * pi) + log(vc0Null)) -
              0.5 / vc0Null * normYtilde2;
            if logLikeNull >= logLikeAlt
              vc0 = vc0Null;
              vc1 = 0;
              logLikeAlt = logLikeNull;
            end
        
            # RLRT test statitic
            statRLRT = 2 * (logLikeAlt - logLikeNull);
        
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statRLRT, evalV, evalAdjV, n, rankX,
                                       WPreSim, device = devices,
                                       nSimPts = nNullSimPts,
                                       nNewtonIter = nNullSimNewtonIter,
                                       test = "eRLRT",
                                       pvalueComputing = pvalueComputings,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
        
            # estimate fixed effects
            if isempty(X)
              b = zeros(0);
            else
              Xrot = UV' * X;
              yrot = UV' * y;
              wt = 1.0 ./ sqrt(vc1 * evalV + vc0);
              Xnew = scale(wt, Xrot);
              ynew = wt .* yrot;
              b = Xnew \ ynew;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          elseif tests == "eScore"
        
            # fit the null model
            b = X \ y;
            r = y - X * b;
            vc0 = norm(r) ^ 2 / n;
            vc1 = 0;
        
            # score test statistic
            #statScore = norm(r' * sqrtV) ^ 2 / norm(r) ^ 2;
            statScore = norm(BLAS.gemv('T', sqrtV, r)) ^ 2 / norm(r) ^ 2;
            #statScore = max(statScore, sum(evalV) / n);
        
            #=
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statScore, evalV, evalAdjV, n, rankX,
                                       WPreSim, test = "eScore",
                                       nSimPts = nNullSimPts,
                                       pvalueComputing = pvalueComputings,
                                       nNewtonIter = nNullSimNewtonIter,
                                       device = devices,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
            =#
        
            if statScore <= sum(evalV) / n
              vc1_pvalue = 1.0;
            else
              w = [evalAdjV - statScore; - statScore];
              dfvec = [ones(Int, rankAdjV); n - rankX - rankAdjV];
              fun(x) = imag(prod((1 - 2 .* w * x * im) .^ (- dfvec / 2)) / x);
              (vc1_pvalue, err) = quadgk(fun, 0, Inf);
              vc1_pvalue = 0.5 + vc1_pvalue / pi;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          end
        end
 #------- 
          (b, vc0List[g], vc1List[g], pvalList[g]) = testvc2();
        end

  elseif ZtZ =="intercept"

        # obtain some orthonormal vectors of space R^n - [QX,QPHI,S]
        # See Golub and Van Loan (1996) Algorithm 12.4.2 on p602

     BLAS.gemm!('T', 'N', 1.0, K, XPhitNullBasis, 0.0, tmpMat); ##zj:: gSNP replaced by K
        ##zj:: the tmpMat here is 
     (UXPhitS, svalXPhitS, VXPhitS) = svd(tmpMat, thin = false);
      # if g == nGrps && grpSize != windowSize
      #   tmpMat = [tmpMat; Array(Float64, windowSize - grpSize,
      #                           size(XPhitNullBasis, 2))];
      # end
     rankXPhitS = countnz(svalXPhitS .> size(XPhitNullBasis, 2) * eps(svalXPhitS[1]));  ##zj:: max(size(XPhitNullBasis, 2), 
    # pXPhitSV = pointer(VXPhitS) + rankXPhitS * size(XPhitNullBasis, 2) * sizeof(Float64);
     pXPhitSV = pointer(VXPhitS)
     BLAS.blascopy!(size(XPhitNullBasis, 2) * rankQZ, pXPhitSV, 1, pSubXPhitSV, 1);
     BLAS.gemm!('N', 'N', 1.0, XPhitNullBasis, subXPhitSV, 0.0, QRes); # QRes---Q2
     # working Y and working non-trivial variance component
     #yWork = (PhiAdjsvd[:U]' * (yShift + KPhiAdj * (QRes' * y))) ./ sqrt(evalPhiAdj);
     BLAS.blascopy!(length(yShift), yShift, 1, tmpvec, 1); # tmpvec---Q1'*y
     BLAS.gemv!('T', 1.0, QRes, y, 0.0, tmpvecQRes);  # tmpvecQRes---Q2'*y
     BLAS.gemv!('N', 1.0, KPhiAdj, tmpvecQRes, 1.0, tmpvec); #tmpvec---Q1'*y+K*Q2'*Y
     #yWork = BLAS.gemv('T', 1.0, PhiAdjsvd[:U], tmpvec) ./ sqrt(evalPhiAdj);
     BLAS.gemv!('T', 1.0, weightedW, tmpvec, 0.0, yWork); 

     BLAS.gemm!('T', 'N', 1.0, QZ, K, 0.0, VWorkSqrt2);  #VWorkSqrt2=Q1'*Kernel
     BLAS.gemm!('T', 'N', 1.0, weightedW, VWorkSqrt2, 0.0, VWorkSqrt);
     Vform="half";
     #scale!(InvSqrtevalPhiAdj, VWorkSqrt);

       #----------------testvc is for eRLRT and eScore(longitudinal)---------------#    

    function testvc(y=yWork,V=VWorkSqrt, X=[],tests=test, devices=device, vcInit = Float64[],
                       bInit = Float64[], WPreSim= [Float64[] Float64[]],
                       pvalueComputings = pvalueComputing,nMMmax=nMMmax)
       
         n = length(y);
         if size(V, 1) != n
           error("vctest:wrongdimV\n", "dimension of V does not match that of X");
         end
       
         # set default maximum MM iteration
         if nMMmax == 0 && tests == "eLRT"
           nMMmax = 10;
         elseif nMMmax == 0 && tests == "eRLRT"
           nMMmax = 1000;
         end
       
         # SVD of X
         if isempty(X)
           rankX = 0;
           X = reshape(X, n, 0);
           X = convert(Array{Float64, 2}, X);
           # LRT is same as RLRT if X is null
           if tests == "eLRT"
             tests = "eRLRT";
           end
         else
           if tests == "eRLRT"
             (UX, svalX) = svd(X, thin = false);
           else
             (UX, svalX) = svd(X);
           end
           rankX = countnz(svalX .> n * eps(svalX[1]));
         end
       
         # eigendecomposition of V
         if Vform == "whole"
           (evalV, UV) = eig(V);
           rankV = countnz(evalV .> n * eps(sqrt(maximum(evalV))));
           sortIdx = sortperm(evalV, rev = true);
           evalV = evalV[sortIdx[1:rankV]];
           UV = UV[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = V;
           end
         elseif Vform == "half"
           (UVfull, tmpevalVfull) = svd(V, thin=false);
           evalVfull = zeros(n);
           pevalVfull = pointer(evalVfull);
           ptmpevalV = pointer(tmpevalVfull);
           BLAS.blascopy!(length(tmpevalVfull), ptmpevalV, 1, pevalVfull, 1);
           rankV = countnz(evalVfull .> n * eps(maximum(evalVfull)));
           evalVfull = evalVfull .^ 2;
           evalV = evalVfull[1:rankV];
           UV = UVfull[:, 1:rankV];
           #if tests == "eLRT" || tests == "eRLRT"
           #  wholeV = *(V, V');
           #end
         elseif Vform == "eigen"
           UV = V.U;
           evalV = V.eval;
           rankV = countnz(V.eval .> n * eps(sqrt(maximum(V.eval))));
           sortIdx = sortperm(V.eval, rev = true);
           V.eval = V.eval[sortIdx[1:rankV]];
           V.U = V.U[:, sortIdx[1:rankV]];
           if tests == "eLRT" || tests == "eRLRT"
             wholeV = *(V.U .* reshape(V.eval, length(V.eval), 1), V.U');
           end
         end
       
         # obtain eigenvalues of (I-PX)V(I-PX)
         if !isempty(X) || tests == "eScore"
           #sqrtV = UV .* sqrt(evalV)';
           sqrtV = similar(UV);
           psqrtV = pointer(sqrtV);
           pUV = pointer(UV);
           BLAS.blascopy!(n*rankV, pUV, 1, psqrtV, 1);
           scale!(sqrtV, sqrt(evalV));
           # scale!(sqrtV, sqrt(abs(evalV)));
         end
         if isempty(X)
           evalAdjV = evalV;
         else
           subUX = Array(Float64, n, rankX);
           psubUX = pointer(subUX);
           pUX = pointer(UX);
           BLAS.blascopy!(n*rankX, pUX, 1, psubUX, 1);
           mat1 = BLAS.gemm('T', 'N', 1.0, subUX, sqrtV);
           mat2 = BLAS.gemm('N', 'N', 1.0, subUX, mat1);
           (UAdjV, evalAdjV) = svd(sqrtV - mat2, thin = false);
           if isempty(evalAdjV)
             evalAdjV = Float64[];
           else
             evalAdjV = evalAdjV[evalAdjV .> n * eps(maximum(evalAdjV))] .^ 2;
           end
         end
         rankAdjV = length(evalAdjV);
       
         srand(1)
         nSimPts = nNullSimPts;
         if size(WPreSim, 1) < rankAdjV
             newSim = randn(rankAdjV - size(WPreSim, 1), nSimPts) .^ 2;
             if isempty(WPreSim)
               WPreSim = newSim;    ##WPreSim=54,10000
             else
               WPreSim = [WPreSim, newSim];
             end
             #windowSize = rankAdjV;
         end
       
         if test == "eRLRT" || test == "eScore"
          windowSize=rankAdjV;
         elseif test == "eLRT"
          windowSize=rankV;
         end
       
         if ZtZ =="intercept"
         nPreRank = min(length(evalPhiAdj), windowSize);
         else
          nPreRank=n-rankX;  #zj::
         end
         # nPreRank= rankAdjV
        
         PrePartialSumW = Array(Float64, nNullSimPts, nPreRank+1);
         PreTotalSumW = Array(Float64, nNullSimPts, nPreRank+1);
         tmpSumVec = Array(Float64, nNullSimPts);
         p1 = pointer(PrePartialSumW);
         BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn), nNullSimPts), 1, p1, 1);
         for i = 1 : nNullSimPts
           #pW = pointer(WPreSim) + (i - 1) * windowSize * sizeof(Float64);
           #tmpSumVec[i] = BLAS.asum(1, pW, 1);
           tmpSumVec[i] = 0.0;
           PreTotalSumW[i, 1] = tmpSumVec[i] + PrePartialSumW[i, 1];
         end
         for j = 1 : nPreRank
           p1 = pointer(PrePartialSumW) + j * nNullSimPts * sizeof(Float64);
           BLAS.blascopy!(nNullSimPts, rand(Distributions.Chisq(tmpn +1 - j), nNullSimPts), 1, p1, 1); ##ZJ:: tmpn-j???
           for i = 1 : nNullSimPts
             tmpSumVec[i] += WPreSim[j, i];
             PreTotalSumW[i, j+1] = tmpSumVec[i] + PrePartialSumW[i, j+1];
           end
         end
       
         if pvalueComputing == "chi2"
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nPtsChi2);
           totalSumW = Array(Float64, nPtsChi2);
           lambda = Array(Float64, 1, nPtsChi2);
           W = Array(Float64, windowSize, nPtsChi2);
           tmpmat0 = Array(Float64, windowSize, nPtsChi2);
           tmpmat1 = Array(Float64, windowSize, nPtsChi2);
           tmpmat2 = Array(Float64, windowSize, nPtsChi2);
           tmpmat3 = Array(Float64, windowSize, nPtsChi2);
           tmpmat4 = Array(Float64, windowSize, nPtsChi2);
           tmpmat5 = Array(Float64, windowSize, nPtsChi2);
           denomvec = Array(Float64, nPtsChi2);
           d1f = Array(Float64, nPtsChi2);
           d2f = Array(Float64, nPtsChi2);
           if test == "eScore"
             simnull = Array(Float64, nPtsChi2);
           else
             simnull = Float64[];
           end
         else
           nPtsChi2 = 300;
           partialSumW = Array(Float64, nNullSimPts);
           totalSumW = Array(Float64, nNullSimPts);
           lambda = Array(Float64, 1, nNullSimPts);
           W = Array(Float64, windowSize, nNullSimPts);
           tmpmat0 = Array(Float64, windowSize, nNullSimPts);
           tmpmat1 = Array(Float64, windowSize, nNullSimPts);
           tmpmat2 = Array(Float64, windowSize, nNullSimPts);
           tmpmat3 = Array(Float64, windowSize, nNullSimPts);
           tmpmat4 = Array(Float64, windowSize, nNullSimPts);
           tmpmat5 = Array(Float64, windowSize, nNullSimPts);
           denomvec = Array(Float64, nNullSimPts);
           d1f = Array(Float64, nNullSimPts);
           d2f = Array(Float64, nNullSimPts);
           if test == "eScore"
             simnull = Array(Float64, nNullSimPts);
           else
             simnull = Float64[];
           end
         end
       #end
       
         # fit the variance component model
         if tests == "eLRT"
       
           # estimates under null model
           bNull = X \ y;
           rNull = y - X * bNull;
           vc0Null = norm(rNull) ^ 2 / n;
           Xrot = UVfull' * X;
           yrot = UVfull' * y;
           loglConst = - 0.5 * n * log(2.0 * pi);
       
           # set starting point
           if isempty(bInit)
             b = copy(bNull);
             r = copy(rNull);
           else
             b = copy(bInit);
             r = y - X * b;
           end
           if isempty(vcInit)
             vc0 = norm(r) ^ 2 / n;
             vc1 = 1;
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
           else
             vc0 = vcInit[1];
             # avoid sticky boundary
             if vc0 == 0
               vc0 = 1e-4;
             end
             vc1 = vcInit[2];
             if vc1 == 0
               vc1 = 1e-4;
             end
             wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
           end
       
           # update residuals according supplied var. components
           r = y - BLAS.gemv('N', X, b);
           rnew = BLAS.gemv('T', UVfull, r);
           #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
           logl = loglConst + sum(log, wt) - 0.5 * sumabs2(rnew .* wt);
       
           nIters = 0;
           #bOld = similar(b);
           #pbOld = pointer(bOld);
           #pb = pointer(b);
           #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
           #denvec = similar(rnew);
           denvec = Float64[];
           for iBlockAscent = 1:nBlockAscent
       
             nIters = iBlockAscent;
       
             # update variance components
             for iMM = 1:nMMmax
               denvec = 1.0 ./ (vc1 * evalVfull + vc0);
               numvec = rnew .* denvec;
               #vc0 = sqrt( (vc0 ^ 2 * sumabs2(numvec) + deltaRes) /
               #            (sumabs(denvec) + (n - rankV) / vc0) );
               vc0 = vc0 * sqrt(sumabs2(numvec) / sumabs(denvec));
               vc1 = vc1 * sqrt( dot(numvec, numvec .* evalVfull) /
                                  sumabs(evalVfull .* denvec) );
               #wt = 1.0 ./ sqrt(vc1 * evalVfull + vc0);
             end
       
             # update fixed effects and residuals
             #pb = pointer(b);
             #BLAS.blascopy!(length(b), pb, 1, pbOld, 1);
             wt = sqrt(denvec);
             Xnew = scale(wt, Xrot);
             ynew = wt .* yrot;
             b = Xnew \ ynew;
             r = y - BLAS.gemv('N', X, b);
             rnew = BLAS.gemv('T', UVfull, r);
             #deltaRes = norm(r) ^ 2 - norm(rnew) ^ 2;
       
             # stopping criterion
             loglOld = logl;
             logl = loglConst + sum(log, wt) -
               0.5 * BLAS.dot(length(evalVfull), rnew .^ 2, 1, denvec, 1);
             if abs(logl - loglOld) < tolX * (abs(logl) + 1.0)
               break
             end
       
           end
       
           # log-likelihood at alternative
           logLikeAlt = logl;
           # log-likelihood at null
           logLikeNull = loglConst - 0.5 * n * log(vc0Null) -
             0.5 * sum(rNull .^ 2) / vc0Null;
           if logLikeNull >= logLikeAlt
             vc0 = vc0Null;
             vc1 = 0;
             logLikeAlt = logLikeNull;
             b = bNull;
             r = rNull;
           end
       
           # LRT test statistic
           statLRT = 2 * (logLikeAlt - logLikeNull);
       
           # obtain p-value for testing vc1=0
           vc1_pvalue = vctestnullsim(statLRT, evalV, evalAdjV, n, rankX,
                                      WPreSim, device = devices,
                                      nSimPts = nNullSimPts,
                                      nNewtonIter = nNullSimNewtonIter,
                                      test = "eLRT",
                                      pvalueComputing = pvalueComputings,
                                      PrePartialSumW = PrePartialSumW,
                                      PreTotalSumW = PreTotalSumW,
                                      partialSumWConst = partialSumWConst,
                                      totalSumWConst = totalSumWConst,
                                      windowSize = windowSize,
                                      partialSumW = partialSumW,
                                      totalSumW = totalSumW,
                                      lambda = lambda, W = W,
                                      nPreRank = nPreRank,
                                      tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                      tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                      tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                      denomvec = denomvec,
                                      d1f = d1f, d2f = d2f, offset = offset,
                                      nPtsChi2 = nPtsChi2, simnull = simnull);
       
           # return values
           return b, vc0, vc1, vc1_pvalue;
       
         elseif tests == "eRLRT"
       
           if isempty(X)
             ytilde = y;
             rankBVB = rankV;
             evalBVB = evalV;
             UBVB = UV;
           else
             # obtain a basis of N(X')
             B = UX[:, rankX+1:end];
             ytilde = B' * y;
       
             # eigen-decomposition of B'VB and transform data
             (UBVB, evalBVB) = svd(B' * sqrtV);
             rankBVB = countnz(evalBVB .> n * eps(maximum(evalBVB)));
             evalBVB = evalBVB[1:rankBVB] .^ 2;
             UBVB = UBVB[:, 1:rankBVB];
           end
           resnew = UBVB' * ytilde;
           normYtilde2 = norm(ytilde) ^ 2;
           deltaRes = normYtilde2 - norm(resnew) ^ 2;
       
           # set initial point
           # TODO: is there better initial point?
           vc0Null = normYtilde2 / length(ytilde);
           if isempty(vcInit)
             vc0 = vc0Null;
             vc1 = 1.0;
           else
             vc0 = vcInit[1];
             vc1 = vcInit[2];
           end
       
      # MM loop for estimating variance components
            nIters = 0;
            #tmpvec = Array(Float64, rankBVB);
            #denvec = Array(Float64, rankBVB);
            #numvec = Array(Float64, rankBVB);
            tmpvec = Float64[];
            for iMM = 1:nMMmax
              nIters = iMM;
              #tmpvec = vc0 + vc1 * evalBVB;
              #tmpvec = evalBVB;
              #numvecSum = 0.0;
              #denvecSum = 0.0;
              #numvecProdSum = 0.0;
              #denvecProdSum = 0.0;
              tmpvec = BLAS.scal(rankBVB, vc1, evalBVB, 1);
              #for i = 1 : rankBVB
              #  tmpvec[i] += vc0;
              #  denvec[i] = 1 / tmpvec[i];
              #  numvec[i] = (resnew[i] * denvec[i]) ^ 2;
              #  numvecSum += numvec[i];
              #  denvecSum += denvec[i];
              #  numvecProdSum += evalBVB[i] * numvec[i];
              #  denvecProdSum += evalBVB[i] * denvec[i];
              #end
              tmpvec += vc0;
              denvec = 1 ./ tmpvec;
              numvec = (resnew .* denvec) .^ 2;
              vcOld = [vc0 vc1];
              vc0 = sqrt( (vc0 ^ 2 * sum(numvec) + deltaRes) /
                           (sum(denvec) + (n - rankX - rankBVB) / vc0) );
              #vc0 = sqrt( (vc0 ^ 2 * numvecSum + deltaRes) /
              #             (denvecSum + (n - rankX - rankBVB) / vc0) );
              vc1 = vc1 * sqrt( sum(evalBVB .* numvec) / sum(evalBVB .* denvec) );
              #vc1 = vc1 * sqrt( numvecProdSum / denvecProdSum );
              # stopping criterion
              if norm([vc0 vc1] - vcOld) <= tolX * (norm(vcOld) + 1)
                break;
              end
            end
        
            # restrictive log-likelihood at alternative
        
            loglConst = - 0.5 * (n - rankX) * log(2 * pi);
            logLikeAlt =  loglConst - 0.5 * sum(log(tmpvec)) -
              0.5 * (n - rankX - rankBVB) * log(vc0) - 0.5 * normYtilde2 / vc0 +
              0.5 * sum(resnew .^ 2 .* (1 / vc0 - 1 ./ (tmpvec)));
            # restrictive log-likelihood at null
            logLikeNull = - 0.5 * (n - rankX) * (log(2 * pi) + log(vc0Null)) -
              0.5 / vc0Null * normYtilde2;
            if logLikeNull >= logLikeAlt
              vc0 = vc0Null;
              vc1 = 0;
              logLikeAlt = logLikeNull;
            end
        
            # RLRT test statitic
            statRLRT = 2 * (logLikeAlt - logLikeNull);
        
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statRLRT, evalV, evalAdjV, n, rankX,
                                       WPreSim, device = devices,
                                       nSimPts = nNullSimPts,
                                       nNewtonIter = nNullSimNewtonIter,
                                       test = "eRLRT",
                                       pvalueComputing = pvalueComputings,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
        
            # estimate fixed effects
            if isempty(X)
              b = zeros(0);
            else
              Xrot = UV' * X;
              yrot = UV' * y;
              wt = 1.0 ./ sqrt(vc1 * evalV + vc0);
              Xnew = scale(wt, Xrot);
              ynew = wt .* yrot;
              b = Xnew \ ynew;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          elseif tests == "eScore"
        
            # fit the null model
            b = X \ y;
            r = y - X * b;
            vc0 = norm(r) ^ 2 / n;
            vc1 = 0;
        
            # score test statistic
            #statScore = norm(r' * sqrtV) ^ 2 / norm(r) ^ 2;
            statScore = norm(BLAS.gemv('T', sqrtV, r)) ^ 2 / norm(r) ^ 2;
            #statScore = max(statScore, sum(evalV) / n);
        
            #=
            # obtain p-value for testing vc1=0
            vc1_pvalue = vctestnullsim(statScore, evalV, evalAdjV, n, rankX,
                                       WPreSim, test = "eScore",
                                       nSimPts = nNullSimPts,
                                       pvalueComputing = pvalueComputings,
                                       nNewtonIter = nNullSimNewtonIter,
                                       device = devices,
                                       PrePartialSumW = PrePartialSumW,
                                       PreTotalSumW = PreTotalSumW,
                                       partialSumWConst = partialSumWConst,
                                       totalSumWConst = totalSumWConst,
                                       windowSize = windowSize,
                                       partialSumW = partialSumW,
                                       totalSumW = totalSumW,
                                       lambda = lambda, W = W,
                                       nPreRank = nPreRank,
                                       tmpmat0 = tmpmat0, tmpmat1 = tmpmat1,
                                       tmpmat2 = tmpmat2, tmpmat3 = tmpmat3,
                                       tmpmat4 = tmpmat4, tmpmat5 = tmpmat5,
                                       denomvec = denomvec,
                                       d1f = d1f, d2f = d2f, offset = offset,
                                       nPtsChi2 = nPtsChi2, simnull = simnull);
            =#
        
            if statScore <= sum(evalV) / n
              vc1_pvalue = 1.0;
            else
              w = [evalAdjV - statScore; - statScore];
              dfvec = [ones(Int, rankAdjV); n - rankX - rankAdjV];
              fun(x) = imag(prod((1 - 2 .* w * x * im) .^ (- dfvec / 2)) / x);
              (vc1_pvalue, err) = quadgk(fun, 0, Inf);
              vc1_pvalue = 0.5 + vc1_pvalue / pi;
            end
        
            # return values
            return b, vc0, vc1, vc1_pvalue;
        
          end
        end
 #------- 
       (b, vc0List[g], vc1List[g], pvalList[g]) = testvc()
     end
      idx +=1;
      results[idx, :] = hcat(vc0List[g], vc1List[g], pvalList[g]);

   end #end for g

#  return results[1:idx, :];

#end

  

    # write output
    fid = open(outFile, "a");
    #for i = 1 : ncores
      for j = 1 : size(results, 1)
        println(fid, results[j, 1], ",", results[j, 2], ",",
                results[j, 3]);
      end
   # end
    close(fid);
 end

end


