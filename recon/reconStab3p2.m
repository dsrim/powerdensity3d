function reconStab3p2(gamma_type,N,const,hmax)
%RECONSTAB3P2
%
% implementation of stabilized 3+2 algorithm

save_fname = getSaveFname('reconStab3p2',gamma_type,const,N); 
save_perp_fname = getSaveFname('reconStab3p2_perp',gamma_type,const,N); 
sol_fname = getSaveFname('AnisoSolsStab3p2',gamma_type,const,N);
if exist(sol_fname,'file')
    disp(['saved file found: ' sol_fname])
else
    computeAnisoSolsStab3p2(gamma_type,const,N,hmax);
end
load(sol_fname)

% compute Gamma on the uniform grid
M = N^3;
x = linspace(-1,1,N);
h = 2/(N-1);
[X,Y,Z]   = ndgrid(x,x,x);
gamma_fun = aniso_conductivity(gamma_type,const);
Gamma     = reshape(gamma_fun(X(:)',Y(:)',Z(:)'),3,3,N^3);

H_fname = getSaveFname('AnisoSolsStab3p2_H',gamma_type,const,N);
if exist(H_fname,'file')
    disp(['saved file found: ' sol_fname])
	load(H_fname)
else

% compute power densities H1, H2, H3, H4
H1 = computeAnisoHtoolbox(Gamma,h,N,...
                          result1, result2, result3, result4,result5);
H2 = computeAnisoHtoolbox(Gamma,h,N,...
                          result21,result2, result3, result5,result6);
H3 = computeAnisoHtoolbox(Gamma,h,N,....
                          result1, result22,result3, result4,result6);
H4 = computeAnisoHtoolbox(Gamma,h,N,....
                          result1, result2, result23,result5,result6);

                      
% compute determinant DU for different triples
detDU1 = computeDetDUToobox(X,Y,Z,  result1,  result2,  result3);
detDU2 = computeDetDUToobox(X,Y,Z, result21,  result2,  result3);
detDU3 = computeDetDUToobox(X,Y,Z,  result1, result22,  result3);
detDU4 = computeDetDUToobox(X,Y,Z,  result1,  result2, result23);
save(H_fname,'H1','H2','H3','H4',...
			 'detDU1','detDU2','detDU3','detDU4');
end

                      
%% Stabilized reconstruction algorithm

% reconstruct normalized \tilde{\gamma}
tMu1 = computeAnisoMuStab(H1);                     
tMu2 = computeAnisoMuStab(H2);                     
tMu3 = computeAnisoMuStab(H3);                     
tMu4 = computeAnisoMuStab(H4);                     

% compute tZ (stabilized version)
[tZ1,tZHOm1] = computeNormalMatStab(tMu1,H1,N);
[tZ2,tZHOm2] = computeNormalMatStab(tMu2,H2,N);
[tZ3,tZHOm3] = computeNormalMatStab(tMu3,H3,N);    
[tZ4,tZHOm4] = computeNormalMatStab(tMu4,H4,N);

save(save_perp_fname,'tZ1','tZ2','tZ3','tZ4',...
                     'tZHOm1','tZHOm2','tZHOm3','tZHOm4')

% compute tB 
[tB1,~,~,condA1] = computePerpBStab(tZ1, tZHOm1, H1);
[tB2,~,~,condA2] = computePerpBStab(tZ2, tZHOm2, H2);
[tB3,~,~,condA3] = computePerpBStab(tZ3, tZHOm3, H3);
[tB4,~,~,condA4] = computePerpBStab(tZ4, tZHOm4, H4);

% compute G parallel to tAs
G1 = computeGStab(tB1,H1);
G2 = computeGStab(tB2,H2);
G3 = computeGStab(tB3,H3);
G4 = computeGStab(tB4,H4);

% compute tgamma from G
detG1 = repmat(reshape(det3(G1),1,1,M),3,3,1);
detG2 = repmat(reshape(det3(G2),1,1,M),3,3,1);
detG3 = repmat(reshape(det3(G3),1,1,M),3,3,1);
detG4 = repmat(reshape(det3(G4),1,1,M),3,3,1);

tgamma1 = G1 ./ (detG1.^(1/3));
tgamma2 = G2 ./ (detG2.^(1/3));
tgamma3 = G3 ./ (detG3.^(1/3));
tgamma4 = G4 ./ (detG4.^(1/3));

% compute error for tgamma
e4tgamma1 = computeAnisopartError(Gamma,tgamma1);
e4tgamma2 = computeAnisopartError(Gamma,tgamma2);
e4tgamma3 = computeAnisopartError(Gamma,tgamma3);
e4tgamma4 = computeAnisopartError(Gamma,tgamma4);

maxe4tgamma1 = max(abs(e4tgamma1(:)));
maxe4tgamma2 = max(abs(e4tgamma2(:)));
maxe4tgamma3 = max(abs(e4tgamma3(:)));
maxe4tgamma4 = max(abs(e4tgamma4(:)));

meane4tgamma1 = mean(abs(e4tgamma1(:)));
meane4tgamma2 = mean(abs(e4tgamma2(:)));
meane4tgamma3 = mean(abs(e4tgamma3(:)));
meane4tgamma4 = mean(abs(e4tgamma4(:)));

disp(['tgamma1 | max error is: ' num2str(maxe4tgamma1) ...
             ' , mean error is: ' num2str(meane4tgamma1) '\n'])
disp(['tgamma2 | max error is: ' num2str(maxe4tgamma2) ...
             ' , mean error is: ' num2str(meane4tgamma2) '\n'])
disp(['tgamma3 | max error is: ' num2str(maxe4tgamma3) ...
             ' , mean error is: ' num2str(meane4tgamma3) '\n'])
disp(['tgamma3 | max error is: ' num2str(maxe4tgamma4) ...
             ' , mean error is: ' num2str(meane4tgamma4) '\n'])


% compute weighted wtgamma
detH1 = repmat(reshape(det3(H1),1,1,M),3,3,1);
detH2 = repmat(reshape(det3(H2),1,1,M),3,3,1);
detH3 = repmat(reshape(det3(H3),1,1,M),3,3,1);
detH4 = repmat(reshape(det3(H4),1,1,M),3,3,1);

w1 = detH1;
w2 = detH2;
w3 = detH3;
w4 = detH4;

wtgamma = w1.*tgamma1 + w2.*tgamma2 + w3.*tgamma3 + w4.*tgamma4;
detfactor = repmat(reshape(det3(wtgamma),1,1,M),3,3,1).^(1/3);
wtgamma = wtgamma ./ detfactor;
e4wtgamma = computeAnisopartError(Gamma,wtgamma);
maxe4wtgamma = max(abs(e4wtgamma(:)));
meane4wtgamma = mean(abs(e4wtgamma(:)));

disp(['wtgamma | max error is: ' num2str(maxe4wtgamma) ...
             ' , mean error is: ' num2str(meane4wtgamma) '\n'])


% reconstruct scaling tau
b1 = computeRSrcStab(G1,H1,tB1);
b2 = computeRSrcStab(G2,H2,tB2);
b3 = computeRSrcStab(G3,H3,tB3);
b4 = computeRSrcStab(G4,H4,tB4);

detH1 = det3(H1); detH1 = repmat(reshape(detH1,1,1,N^3),3,3,1);
detH2 = det3(H2); detH2 = repmat(reshape(detH2,1,1,N^3),3,3,1);
detH3 = det3(H3); detH3 = repmat(reshape(detH3,1,1,N^3),3,3,1);
detH4 = det3(H4); detH4 = repmat(reshape(detH4,1,1,N^3),3,3,1);

GdetH1 = G1.*detH1; 
GdetH2 = G2.*detH2; 
GdetH3 = G3.*detH3; 
GdetH4 = G4.*detH4; 

b = cramersrule3(GdetH1 + GdetH2 + GdetH3 + GdetH4,(b1 + b2 + b3 + b4)')';

divb = computeDivVec(b,h);

% compute true scaling, Tau
trueTau = det3(Gamma).^(1/3);

% set BC as log(tau) for the Poisson problem on the uniform grid
[xbd0,xbd1,ybd0,ybd1,zbd0,zbd1] = getBdry(N);   
logtau = zeros(N,N,N);
logtau(xbd0) = log(trueTau(xbd0)); 
logtau(xbd1) = log(trueTau(xbd1));
logtau(ybd0) = log(trueTau(ybd0));    
logtau(ybd1) = log(trueTau(ybd1));
logtau(zbd0) = log(trueTau(zbd0));    
logtau(zbd1) = log(trueTau(zbd1));
 
% solve Poisson problem with RHS and BC
logtau = solvePoisson(logtau(:),N,divb(:));
logtau = reshape(logtau,N,N,N);
tau = exp(logtau);

% compute errors
e4tau = trueTau(:) - tau(:);
maxe4tau = max(abs(e4tau(:)));
meane4tau = mean(abs(e4tau(:)));
display(['tau  | max error is: ' num2str(maxe4tau) ...
            ' , mean error is: ' num2str(meane4tau) ])

% save outputs
save(save_fname, 'G1', 'G2', 'G3', 'G4',...
			'tau', 'Gamma', ...
			'H1', 'H2', 'H3', 'H4', ...
            'tgamma1','e4tgamma1',...
            'tgamma2','e4tgamma2',...
            'tgamma3','e4tgamma3',...
            'tgamma4','e4tgamma4',...
            'e4tau',...
            'detDU1','detDU2','detDU3','detDU4',...
			'condA1','condA2','condA3','condA4');


end

