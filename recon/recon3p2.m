function recon3p2(gamma_type,N,const,hmax)
%RECON3P2
%
% implementation of 3+2 algorithm

save_output = true;

sol_fname = getSaveFname('AnisoSols3p2',gamma_type,const,N);
if exist(sol_fname,'file')
    disp(['saved file found: ' sol_fname])
else
    computeAnisoSols3p2(gamma_type,const,N,hmax);
end
load(sol_fname)

x = linspace(-1,1,N);
h = 2/(N-1);
[X,Y,Z]   = ndgrid(x,x,x);
gamma_fun = aniso_conductivity(gamma_type,const);
Gamma     = reshape(gamma_fun(X(:)',Y(:)',Z(:)'),3,3,N^3);

% compute power densities H1, H2, H3, H4
H = computeAnisoHtoolbox(Gamma,h,N,...
                         result1, result2, result3, result4, result5);
                      
% compute determinant of DU 
detDU = computeDetDUToobox(X,Y,Z, result1, result2, result3);
   
Mu = computeAnisoMu(H);

[Z,ZHOm] = computeNormalMat(Mu,H,N);
[B, detvalH,detvalB,condA] = computePerpB(Z, ZHOm, H);

tgamma = computeTgamma(B,H);        % compute tgamma
e4tgamma = computeAnisopartError(Gamma,tgamma);
maxe4tgamma = max(abs(e4tgamma(:)));
meane4tgamma = mean(abs(e4tgamma(:)));
display(['tgamma | max error is: ' num2str(maxe4tgamma) ...
         ' , mean error is: ' num2str(meane4tgamma) ])

%   reconstruct \tau
b = computeRSrc(tgamma,H,B);
b = cramersrule3(tgamma,b')';
divb = computeDivVec(b,h);

% compute true scaling: Tau
trueTau = det3(Gamma).^(1/3);

% set BC in log version for the Poisson problem
[xbd0,xbd1,ybd0,ybd1,zbd0,zbd1] = getBdry(N); 
logtau = zeros(N,N,N);
logtau(xbd0) = log(trueTau(xbd0)); 
logtau(xbd1) = log(trueTau(xbd1));
logtau(ybd0) = log(trueTau(ybd0));    
logtau(ybd1) = log(trueTau(ybd1));
logtau(zbd0) = log(trueTau(zbd0));    
logtau(zbd1) = log(trueTau(zbd1));
 
% solve Poisson problem with the computed RHS
logtau = solvePoisson(logtau(:),N,divb(:)); 
logtau = reshape(logtau,N,N,N);
tau = exp(logtau);

% compute errors
% e4tgamma : error for anisotropic part \tilde{\gamma}
% e4gamma : error for the scaling part \tau
[e4tgamma,e4gamma,e4tau] = computeAnisoError(Gamma,tgamma,tau);
maxe4tau = max(abs(e4tau(:)));
meane4tau = mean(abs(e4tau(:)));
display(['tau | max  error is: ' num2str(maxe4tau) ...
            ' , mean error is: ' num2str(meane4tau) ])

save_fname = getSaveFname('recon3p2',gamma_type,const,N);

% save outputs
if save_output
    save(save_fname, 'B', 'tau', 'Gamma', 'H','detDU','tgamma','tau',...
                     'e4tgamma','e4gamma','e4tau');
end

end

