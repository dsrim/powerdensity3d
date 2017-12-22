
% compute weighted approximation to tgamma
N = 64; M = N^3;

% weight = 'detH';
weight = 'norm9';
% weight = 'indicateH';
% weight = 'indicateCondA';

%%
switch (weight)
    case 'detH'
        detH1 = repmat(reshape(det3(H1),1,1,M),3,3,1);
        detH2 = repmat(reshape(det3(H2),1,1,M),3,3,1);
        detH3 = repmat(reshape(det3(H3),1,1,M),3,3,1);
        detH4 = repmat(reshape(det3(H4),1,1,M),3,3,1);
        
        w1 = detH1;
        w2 = detH2;
        w3 = detH3;
        w4 = detH4;
        
    case 'indicateH'
        detH1 = repmat(reshape(det3(H1),1,1,M),3,3,1);
        detH2 = repmat(reshape(det3(H2),1,1,M),3,3,1);
        detH3 = repmat(reshape(det3(H3),1,1,M),3,3,1);
        detH4 = repmat(reshape(det3(H4),1,1,M),3,3,1);
        
        w1 = detH1;
        w2 = detH2;
        w3 = detH3;
        w4 = detH4;
        
        [v,ii] = max([w1(:),w2(:),w3(:),w4(:)]');
        ii = ii';
        
        w1 = false(size(w1));
        w2 = false(size(w2));
        w3 = false(size(w3));
        w4 = false(size(w4));
        
        w1(ii == 1) = true;
        w2(ii == 2) = true;
        w3(ii == 3) = true;
        w4(ii == 4) = true;       
        

    case 'indicateCondA'
        ii = condA1 < condA2;
        
        w1 = repmat(reshape(ii,1,1,M),3,3,1);
        w2 = repmat(reshape(~ii,1,1,M),3,3,1);
        
    case 'norm9'

        w1 = 1./repmat(reshape(norm9(tgamma1),1,1,M),3,3,1);
        w2 = 1./repmat(reshape(norm9(tgamma2),1,1,M),3,3,1);
        w3 = 1./repmat(reshape(norm9(tgamma3),1,1,M),3,3,1);
        w4 = 1./repmat(reshape(norm9(tgamma4),1,1,M),3,3,1);
        
%         [v,ii] = max([w1(:),w2(:),w3(:),w4(:)]');
%         ii = ii';
%         
%         w1 = false(size(w1));
%         w2 = false(size(w2));
%         w3 = false(size(w3));
%         w4 = false(size(w4));
%         
%         w1(ii == 1) = true;
%         w2(ii == 2) = true;
%         w3(ii == 3) = true;
%         w4(ii == 4) = true;

           
        
end




wtgamma = w1.*tgamma1 + w2.*tgamma2 + w3.*tgamma3 + w4.*tgamma4;
detfactor = repmat(reshape(det3(wtgamma),1,1,M),3,3,1).^(1/3);
wtgamma = wtgamma ./ detfactor;

% - compute true scaling: Tau
trueTau = det3(Gamma).^(1/3);
trueTgamma = Gamma ./ repmat(reshape(trueTau,1,1,M),3,3,1);


[e4wtgamma, e4gamma, e4tau] = computeAnisoError(Gamma,wtgamma,tau);
[~       ,e4gamma1,e4tau1] = computeAnisoError(Gamma,tgamma1,tau);

%%

 figure(1)  ; plotTensor3D(wtgamma,1);

%%
figure; 
nslices = 1;
norme1 = reshape(log10(norm9(e4wtgamma)),N,N,N);
norme2 = reshape(log10(norm9(e4tgamma2)),N,N,N);
vmax = max([ norme1(:); norme2(:)]);
vmin = min([ norme1(:); norme2(:)]);
subplot(131); ezplot3D(norme1,nslices); alpha('color');
title(weight);
subplot(132); ezplot3D(norme2,nslices); alpha('color');
 title('e4tgamma1');
subplot(133); ezplot3D(norme2 - norme1,1); %alpha('color');
title('difference')
set(gcf,'Position',[0,400,1000,300])

%%
figure; 
nslices = 1;
norme1 = reshape(norm9(e4wtgamma)./norm9(e4wtgamma + wtgamma),N,N,N);
norme2 = reshape(norm9(e4tgamma1)./norm9(e4wtgamma + wtgamma),N,N,N);
vmax = max([ norme1(:); norme2(:)]);
vmin = min([ norme1(:); norme2(:)]);
subplot(131); ezplot3D(norme1,nslices); alpha('color');
title(weight);
subplot(132); ezplot3D(norme2,nslices); alpha('color');
 title('e4tgamma1');
subplot(133); ezplot3D(norme2 - norme1,1); %alpha('color');
title('difference')
set(gcf,'Position',[0,400,1000,300])

%%
figure;
subplot(131)
ezplot3D(reshape(tau.*(abs(tau) < 12),N,N,N),3); alpha('color')
% caxis([min(tau(:)),max(tau(:))])
title('approx tau')
subplot(132)
ezplot3D(reshape(trueTau,N,N,N),3); alpha('color')
set(gcf,'Position',[0,200,300,300])
% caxis([min(tau(:)),max(tau(:))])
title('true tau')
subplot(133)
ezplot3D(reshape(abs(e4tau),N,N,N),1); alpha('color')
set(gcf,'Position',[0,200,300,300])
title('e4tau')


figure;
ezplot3D(reshape(abs(e4tau),N,N,N),1); alpha('color')
set(gcf,'Position',[0,200,300,300])
title('e4tau')
figure;
ezplot3D(reshape(norm9(e4gamma)./norm9(Gamma),N,N,N),3); alpha('color')
set(gcf,'Position',[0,200,300,300])

%%
figure; 
nslices = 3;
normTgamma = reshape(norm9(tgamma1),N,N,N);
ezplot3D(normTgamma,nslices); alpha('color'); title(weight);
set(gcf,'Position',[0,400,500,500])

%%

h = 2/(N-1);

tMu1 = computeAnisoMuStab(H1);                 
tMu2 = computeAnisoMuStab(H2);                 
tMu3 = computeAnisoMuStab(H3);                 


[tZ1,tZHOm1] = computeNormalMatStab(tMu1,H1,N);    
[tZ2,tZHOm2] = computeNormalMatStab(tMu2,H2,N);    
[tZ3,tZHOm3] = computeNormalMatStab(tMu3,H3,N);  

[tB1,detvalH,detvalB1,condA1] = computePerpBStab(tZ1, tZHOm1, H1);
[tB2,      ~,detvalB2,condA2] = computePerpBStab(tZ2, tZHOm2, H2);
[tB3,      ~,detvalB3,condA3] = computePerpBStab(tZ3, tZHOm3, H3);

G1 = computeGStab(tB1,H1);        
G2 = computeGStab(tB2,H2); 
G3 = computeGStab(tB3,H3);


% - reconstruct \tau
b1 = computeRSrcStab(G1,H1,tB1);
b2 = computeRSrcStab(G2,H2,tB2);
b3 = computeRSrcStab(G3,H3,tB3);

detH1 = det3(H1); detH1 = repmat(reshape(detH1,1,1,N^3),3,3,1);
detH2 = det3(H2); detH2 = repmat(reshape(detH2,1,1,N^3),3,3,1);
detH3 = det3(H3); detH3 = repmat(reshape(detH3,1,1,N^3),3,3,1);

GdetH1 = G1.*detH1; 
GdetH2 = G2.*detH2; 
GdetH3 = G3.*detH3; 

b = cramersrule3(GdetH1 + GdetH2 + GdetH3,(b1 + b2 + b3)')';

divb = computeDivVec(b,h);
%%
% - compute true scaling: Tau
trueTau = det3(Gamma).^(1/3);


%%
% - set BC in log version for the Poisson problem
[xbd0,xbd1,ybd0,ybd1,zbd0,zbd1] = getBdry(N);   % helper function that

% - set BC in log version for the Poisson problem
logtau = zeros(N,N,N);
logtau(xbd0) = log(trueTau(xbd0)); 
logtau(xbd1) = log(trueTau(xbd1));
logtau(ybd0) = log(trueTau(ybd0));    
logtau(ybd1) = log(trueTau(ybd1));
logtau(zbd0) = log(trueTau(zbd0));    
logtau(zbd1) = log(trueTau(zbd1));
 
% - solve anisotropic Poisson problem
Id = getGamma(N,'id');

logtau = solveAnisoPoisson(logtau(:),Id,N,divb(:));
logtau = reshape(logtau,N,N,N);
tau = exp(logtau);

% -- Compute errors
% e4tgamma : error for anisotropic part \tilde{\gamma}
% e4gamma : error for the scaling part \tau

% - compute errors

e4tau = trueTau(:) - tau(:);
maxe4tau = max(abs(e4tau(:)));
meane4tau = mean(abs(e4tau(:)));
display(['tau  | max error is: ' num2str(maxe4tau) ...
         ' , mean error is: ' num2str(meane4tau) ])


