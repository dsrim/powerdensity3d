clear;

gamma_name = 'anisotori';
const = 2;
N = 16;

save_fname = getSaveFname('reconStab3p2',gamma_name,const,N);
load(save_fname);

% set plot size
img_height = 325;
img_width = 300;

detDU = reshape(detDU,N,N,N); %detDU = permute(detDU,[2 1 3]);

f1 = figure(1); 
ezplot3D(detDU,1);
cmap0 = flipud(colormap(parula)); colormap(cmap0)
alpha(0.7)
adjustFigures(f1)
print(f1,'detDU.png','-dpng','-r150');



%%
trueTau = tau(:) + e4tau;
f1 = figure(1); clf; 
ezplot3D(reshape(tau,N,N,N),1); alpha(0.7)
caxis([min(trueTau(:)),max(trueTau(:))])
colormap(cmap0)
adjustFigures(f1)
print(f1,'tau2.png','-dpng','-r150');

f1 = figure(1); clf; 

ezplot3D(reshape(trueTau,N,N,N),1); alpha(0.7)
caxis([min(trueTau(:)),max(trueTau(:))])
colormap(cmap0)
adjustFigures(f1)
print(f1,'trueTau2.png','-dpng','-r150');

f1 = figure(1); clf; 
ezplot3D(log10(reshape(abs(e4tau(:))./abs(trueTau(:)),N,N,N)),1); alpha(0.7)
colormap(cmap0)
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'rele4tau2.png','-dpng','-r150');

%%

f1 = figure(1); clf; 
norme2 = reshape(log10(norm9(e4tgamma)),N,N,N);
ezplot3D(norme2,1); alpha(0.7);
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'e4tgamma.png','-dpng','-r150');

%%
[~, e4gamma, ~] = computeAnisoError(Gamma,tgamma,tau);

f1 = figure(1); clf;
ezplot3D(reshape(log10(norm9(e4gamma)./norm9(Gamma)),N,N,N),1); alpha(0.7)
colorbarLogscale(f1);
adjustFigures(f1)
print(f1,'rele4gamma2.png','-dpng','-r150');


%%
norm9e4tgamma = norm9(e4tgamma);
norm9tgamma = norm9(e4tgamma + tgamma);
rele4tgamma = norm9e4tgamma./norm9tgamma;

display(['max eror = ' num2str(max(abs(norm9e4tgamma(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(norm9e4tgamma(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(norm9e4tgamma(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(rele4tgamma),'%10.8f')])
display(['avg rel error = ' num2str(mean(rele4tgamma),'%10.8f')])

%%

norm9e4gamma = norm9(e4gamma);
norm9gamma = norm9(Gamma);
rele4gamma = norm9e4gamma./norm9gamma;

disp('-')
display(['max eror = ' num2str(max(abs(norm9e4gamma(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(norm9e4gamma(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(norm9e4gamma(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(rele4gamma),'%10.8f')])
display(['avg rel error = ' num2str(mean(rele4gamma),'%10.8f')])


%%


disp('-')
display(['max eror = ' num2str(max(abs(e4tau(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(e4tau(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(e4tau(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(abs(e4tau(:))./abs(tau(:))),'%10.8f')])
display(['avg rel error = ' num2str(mean(abs(e4tau(:))./abs(tau(:))),'%10.8f')])

