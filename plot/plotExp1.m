addpath(genpath('../'));

sigma_name = 'gaussians';
const = 1;
N = 64;
nslices = 3;

fname = getSaveFname('reconSigma',sigma_name,const,N);
load(fname);

% set plot size and colormap
img_height = 450;
img_width = 400;
cmap0 = colormap('parula');

%%
e4sigma = reshape(e4sigma,N,N,N);
trueSigma = sigma + e4sigma;
f1 = figure(1);
clf;
hold on;
ezplot3D(trueSigma,nslices,'meshgrid')
view(3);
axis equal
title('')



%% 
colormap(flipud(cmap0)); alpha(0.5);
f1.Children(1).Position(1) = 0.9;
f1.Children(1).Position(4) = 0.7;
xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
adjustFigures(f1)
print('scalar_tori_slices','-dpng','-r150');


%%
f1 = figure(1);
clf
w0 = -1; w1 = 1;
trueSigma = sigma + e4sigma;
x = linspace(w0,w1,N);
[X,Y,Z] = meshgrid(x,x,x);  % slice() accepts only meshgrids

p = patch(isosurface(X, Y, Z, trueSigma, 2.0));
isonormals(X,Y,Z,sigma, p)
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
p.FaceColor = 'blue';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; 
lighting phong

axis equal
print('scalar_tori_isosurface','-dpng','-r150');


%%
sol_fname = getSaveFname('IsoSols',sigma_name,const,N);
load(sol_fname)

h3 = figure(1); 
clf;
vmin = min(detDU(:)); vmax = max(detDU(:));
hold on; 
detDU = reshape(detDU,N,N,N);
ezplot3D(detDU,nslices,'meshgrid'); alpha(0.5)

caxis([vmin,vmax]);
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
colorbar()
axis equal
title('')
view(3)

%%
h3.Children(1).Position(1) = 0.9;
h3.Children(1).Position(4) = 0.8;
colormap(cmap0)
% set(h3,'Position',[100,300,img_height,img_width])
% saveas(gcf,'scalar_tori.png')
% h.PaperPositionMode = 'auto';
adjustFigures(h3);
print('scalar_detDU','-dpng','-r150');
% ezplot3D(detDU,3);
% cmap = colormap(hot)
% alpha('scaled')

%%

h4 = figure(4); 
clf;
vmin = min(trueSigma(:)); vmax = max(trueSigma(:));
hold on; 

ezplot3D(sigma,nslices,'meshgrid'); 
colormap(flipud(cmap0)); alpha(0.8)
% contourslice(X,Y,Z,detDU,0,0,0,linspace(vmin,vmax,16))
caxis([vmin,vmax]);

% xlabel('x'); ylabel('y'); zlabel('z');
% xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);

axis equal
title('')
view(3)


%%
% h4.Children(1).Position(1) = 0.9;
% h4.Children(1).Position(4) = 0.8;


% set(h4,'Position',[100,300,img_height,img_width])
% saveas(gcf,'scalar_tori.png')
% h.PaperPositionMode = 'auto';
adjustFigures(h4);
print('sigma1','-dpng','-r150');
% ezplot3D(detDU,3);


%%

re4sigma = e4sigma./trueSigma;
h4 = figure(4); close(h4); h4 = figure(4);
clf;
% vmin = min(re4sigma(:)); vmax = max(re4sigma(:));
hold on; 
ezplot3D(log10(abs(re4sigma)),nslices,'meshgrid'); alpha(0.5)
% contourslice(X,Y,Z,detDU,0,0,0,linspace(vmin,vmax,16))
% caxis([vmin,vmax]);
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
colorbar()
axis equal
title('')
view(3)


%%

h4.Children(1).Position(1) = 0.9;
h4.Children(1).Position(4) = 0.8;

colormap(flipud(cmap0)); alpha(0.5);
set(h4,'Position',[100,300,img_height,img_width])
colorbarLogscale(h4)
% saveas(gcf,'scalar_tori.png')
% h.PaperPositionMode = 'auto';
adjustFigures(h4);
print('sigma_rele1','-dpng','-r150');
% ezplot3D(detDU,3);
% cmap = colormap(hot)
% alpha('scaled')

%%

disp('-')
display(['max eror = ' num2str(max(abs(e4sigma(:))),'%10.8f')])
display(['L2 error = ' num2str(norm(e4sigma(:))*(2/(N-1))^(3/2),'%10.8f')])
display(['L1 error = ' num2str(norm(e4sigma(:),1)*(2/(N-1))^(3),'%10.8f')])
display(['max rel error = ' num2str(max(abs(e4sigma(:))./abs(trueSigma(:))),'%10.8f')])
display(['avg rel error = ' num2str(mean(abs(e4sigma(:))./abs(trueSigma(:))),'%10.8f')])


%%

img_height = 325;
img_width = 300;

if 0
f5 = figure(5);
clf;
xlim([-1,1])
ylim([-1,1])

axis equal
 
n0 = 4;


for k = 0:3
    clf;
    if k == 0
         components = [1,2,3];
    else
         components = [k];
    end
    for i = 64:64
        for j = 1:n0:N
            display(['j = ' num2str(j)])

            pos = [i,j,32];
            plotEV(pos,components,gamma_name,N,n0,'true')
        end
    end
    
    if k == 0
       zlim([-1,1])
       ylim([-.1,.1])
    end
    title('')
    set(f5,'Position',[0,0,img_width,img_height]);
    f5.Children(1).Position(1) = 0.05;
    f5.Children(1).Position(2) = 0.05;
    f5.Children(1).Position(3) = 0.95;
    f5.Children(1).Position(4) = 0.95;
    print(['Rtrue_component_' num2str(k)],'-dpng','-r400');
end

%%

f5 = figure(5);
clf;
xlim([-1,1])
ylim([-1,1])

axis equal
 
N = 128;
n0 = 4;

for k = 0:3
    clf;
    if k == 0
         components = [1,2,3];
    else
         components = [k];
    end
    for i = 64:64
        for j = 1:n0:N
            display(['j = ' num2str(j)])

            pos = [i,j,32];
            plotEV(pos,components,gamma_name,N,n0,'recon')
        end
    end
    
    if k == 0
       zlim([-1,1])
       ylim([-.1,.1])
    end
    title('')
    
    set(f5,'Position',[0,0,img_width,img_height]);
    f5.Children(1).Position(1) = 0.05;
    f5.Children(1).Position(2) = 0.05;
    f5.Children(1).Position(3) = 0.95;
    f5.Children(1).Position(4) = 0.95;
    print(['R_component_' num2str(k)],'-dpng','-r400');
end

end
