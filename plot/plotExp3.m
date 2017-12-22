clear;

N = 128;
gamma_type = 'aniso_tori';
const = 20;

getSaveFname('reconStab3p2',gamma_type,const,N);

%%
M = N^3;
xx = linspace(-1,1,N);
[X,Y,Z] = meshgrid(xx,xx,xx);
cmap = flipud(colormap(parula));

%%
detDU = detDU1;
detDU = reshape(detDU,N,N,N); detDU = permute(detDU,[2 1 3]);

f1 = figure(1); 
ezplot3D(detDU,1);
colormap(cmap)
adjustFigures(f1)
print(f1,'detDU1.png','-dpng','-r150');

%%
f1 = figure(1); clf;
minmin = min(detDU(:));
contourslice(X,Y,Z,detDU,0,0,0,(-1:0.05:0.0)*abs(minmin))
xlabel('x'); ylabel('y'); zlabel('z'); colorbar
shading flat; axis equal; axis([-1 1 -1 1 -1 1]); 
adjustFigures(f1)
print(f1,'detDU1_cslice.png','-dpng','-r150');

f1 = figure(1); clf;
detDU = detDU2;
detDU = reshape(detDU,N,N,N); detDU = permute(detDU,[2 1 3]);
ezplot3D(detDU,1);
caxis([min(detDU1(:)), max(detDU1(:))]);
adjustFigures(f1)
print(f1,'detDU2.png','-dpng','-r150');

f1 = figure(1); clf;
minmin = min(min(min(detDU))); detDU = permute(detDU,[2 1 3]);
contourslice(X,Y,Z,detDU,0,0,0,(-1:0.025:0.0)*abs(minmin))
xlabel('x'); ylabel('y'); zlabel('z'); colorbar
shading flat; axis equal; axis([-1 1 -1 1 -1 1]); 
adjustFigures(f1)
print(f1,'detDU2_cslice.png','-dpng','-r150');


f1 = figure(1); clf;
detDU = detDU3;
detDU = reshape(detDU,N,N,N); detDU = permute(detDU,[2 1 3]);
ezplot3D(detDU,1);
caxis([min(detDU1(:)), max(detDU1(:))]);
adjustFigures(f1)
print(f1,'detDU3.png','-dpng','-r150');

f1 = figure(1); clf;
minmin = min(min(min(detDU))); detDU = permute(detDU,[2 1 3]);
contourslice(X,Y,Z,detDU,0,0,0,(-1:0.025:0.0)*abs(minmin))
xlabel('x'); ylabel('y'); zlabel('z'); colorbar
shading flat; axis equal; axis([-1 1 -1 1 -1 1]);
adjustFigures(f1)
print(f1,'detDU3_cslice.png','-dpng','-r150');


f1 = figure(1); clf;
detDU = detDU4;
detDU = reshape(detDU,N,N,N); detDU = permute(detDU,[2 1 3]);
ezplot3D(detDU,1);
caxis([min(detDU1(:)), max(detDU1(:))]);
shading flat; axis equal; axis([-1 1 -1 1 -1 1]); 
adjustFigures(f1)
print(f1,'detDU4.png','-dpng','-r150');

f1 = figure(1); clf;
minmin = min(min(min(detDU))); detDU = permute(detDU,[2 1 3]);
contourslice(X,Y,Z,detDU,0,0,0,(-1:0.025:0.0)*abs(minmin))
xlabel('x'); ylabel('y'); zlabel('z'); colorbar
shading flat; axis equal; axis([-1 1 -1 1 -1 1]); 
adjustFigures(f1)
print(f1,'detDU4_cslice.png','-dpng','-r150');

%%

f1 = figure(1); clf;
xx = linspace(-1,1,N);
[X,Y,Z] = meshgrid(xx,xx,xx);
detDU = detDU1.^2 + detDU2.^2 + detDU3.^2 + detDU4.^2;
detDU = reshape(detDU,N,N,N); detDU = permute(detDU,[2 1 3]);
ezplot3D(detDU,1); alpha(0.7)
adjustFigures(f1)
print(f1,'detDUsqsum.png','-dpng','-r150');

% subplot(122);
% minmin = min(min(min(detDU)));

f1 = figure(1); clf; detDU = permute(detDU,[2 1 3]);
contourslice(X,Y,Z,detDU,0,0,0,(0:0.025:1.))
xlabel('x'); ylabel('y'); zlabel('z'); colorbar
%alpha(0.7);
shading flat; axis equal; axis([-1 1 -1 1 -1 1]); 
adjustFigures(f1)
print(f1,'detDUsqsum_cslice.png','-dpng','-r150');
% title(['|detDU1|^2 + |detDU2|^2 + |detDU3|^2 + |detDU4|^2, min =' num2str(minmin)])
% min = 0.3432




%% compute weighted approximation to tgamma

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
        
          
end

wtgamma = w1.*tgamma1 + w2.*tgamma2 + w3.*tgamma3 + w4.*tgamma4;
detfactor = repmat(reshape(det3(wtgamma),1,1,M),3,3,1).^(1/3);
wtgamma = wtgamma ./ detfactor;

% compute true scaling: Tau
trueTau = det3(Gamma).^(1/3);
trueTgamma = Gamma ./ repmat(reshape(trueTau,1,1,M),3,3,1);


[e4wtgamma, e4gamma,  e4tau] = computeAnisoError(Gamma,wtgamma,tau);
[~       , e4gamma1, e4tau1] = computeAnisoError(Gamma,tgamma1,tau);

%%

nslices = 1;
norme1 = reshape(log10(norm9(e4wtgamma)),N,N,N);
% vmax = max([ norme1(:); norme2(:)]);
% vmin = min([ norme1(:); norme2(:)]);
f1 = figure(1); clf; 
ezplot3D(norme1,nslices); 
colormap(cmap)
alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,['e4' weight '.png'],'-dpng','-r150');

f1 = figure(1); clf; 
norme2 = reshape(log10(norm9(e4tgamma1)),N,N,N);
ezplot3D(norme2,nslices); alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'e4tgamma1.png','-dpng','-r150');

f1 = figure(1); clf; 
norme2 = reshape(log10(norm9(e4tgamma2)),N,N,N);
ezplot3D(norme2,nslices); alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'e4tgamma2.png','-dpng','-r150');

f1 = figure(1); clf; 
norme2 = reshape(log10(norm9(e4tgamma3)),N,N,N);
ezplot3D(norme2,nslices); alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'e4tgamma3.png','-dpng','-r150');

f1 = figure(1); clf; 
norme2 = reshape(log10(norm9(e4tgamma4)),N,N,N);
ezplot3D(norme2,nslices); alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'e4tgamma4.png','-dpng','-r150');




%%
f1 = figure(1); clf; 
ezplot3D(reshape(tau,N,N,N),1); alpha(0.7)
caxis([min(trueTau(:)),max(trueTau(:))])
colormap(cmap)
adjustFigures(f1)
print(f1,'tau3.png','-dpng','-r150');

f1 = figure(1); clf; 
ezplot3D(reshape(trueTau,N,N,N),1); alpha(0.7)
caxis([min(trueTau(:)),max(trueTau(:))])
colormap(cmap)
adjustFigures(f1)
print(f1,'trueTau3.png','-dpng','-r150');


f1 = figure(1); clf; 
ezplot3D(reshape(log10(abs(e4tau)./abs(trueTau')),N,N,N),1); alpha(0.7)
colormap(cmap)
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'rele4tau3.png','-dpng','-r150');

f1 = figure(1); clf;
ezplot3D(reshape(log10(norm9(e4gamma)./norm9(Gamma)),N,N,N),1); alpha(0.7)
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,['rele4gamma3' weight '.png'],'-dpng','-r150');

% %%
% figure; 
% nslices = 3;
% normTgamma = reshape(norm9(tgamma1),N,N,N);
% ezplot3D(normTgamma,nslices); alpha(0.7); title(weight);
% set(gcf,'Position',[0,400,img_height,500])

%%

disp('- tau')
display(['max eror = ' num2str(max(abs(e4tau(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(e4tau(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(e4tau(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(abs(e4tau(:))./abs(tau(:))),'%10.8f')])
display(['avg rel error = ' num2str(mean(abs(e4tau(:))./abs(tau(:))),'%10.8f')])

%%
norm9e4wtgamma = norm9(e4wtgamma);
norm9wtgamma   = norm9(e4wtgamma + wtgamma);
rele4wtgamma   = norm9e4wtgamma./norm9wtgamma;

display(['- tgamma, weight = ' weight])
display(['max eror = ' num2str(max(abs(norm9e4wtgamma(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(norm9e4wtgamma(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(norm9e4wtgamma(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(rele4wtgamma),'%10.8f')])
display(['avg rel error = ' num2str(mean(rele4wtgamma),'%10.8f')])

%%

norm9e4gamma = norm9(e4gamma);
norm9gamma   = norm9(Gamma);
rele4gamma   = norm9e4gamma./norm9gamma;

disp('- gamma')
display(['max eror = ' num2str(max(abs(norm9e4gamma(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(norm9e4gamma(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(norm9e4gamma(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(rele4gamma),'%10.8f')])
display(['avg rel error = ' num2str(mean(rele4gamma),'%10.8f')])
