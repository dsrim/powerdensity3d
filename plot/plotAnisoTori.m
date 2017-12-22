
addpath(genpath('../'))
run('../../MVIRT/initMVIRT.m')


a = 280;
img_width = a / 300*325;
img_height = a;

%% Initialization
cmap=flipud(colormap(parula));

N = 128;
g2 = myC('anisotori',2);
xx = linspace(-1,1,N);
[X,Y,Z] = ndgrid(xx,xx,xx);
Gamma2 = g2(X(:)',Y(:)',Z(:)');
Gamma2 = reshape(Gamma2,3,3,N^3);

tau = det3(Gamma2);
tau = reshape(tau,N,N,N);


npts = 40;    % number of ellipsoid points


%%

f1 = figure(1);
ezplot3D(tau,1);
colormap(cmap)
alpha(0.7)
adjustFigures(f1);
print(gcf,'tau_slices1.png','-dpng','-r150');

%%
N = 40;
g2 = myC('anisotori',2);
xx = linspace(-1,1,N);
[X,Y,Z] = ndgrid(xx,xx,xx);
Gamma2 = g2(X(:)',Y(:)',Z(:)');
Gamma2 = reshape(Gamma2,3,3,N^3);

Gamma_slice1 = Gamma2(:,:,(X == xx(floor(N/2))));
detGamma_slice1 = det3(Gamma_slice1);
cvals1 = getColorSPD(detGamma_slice1,N,cmap);
Gamma_slice1 = reshape(Gamma_slice1,3,3,1,N,N);
cvals1 = reshape(cvals1,1,N,N,3);

Gamma_slice2 = Gamma2(:,:,(Y == xx(floor(N/2))));
detGamma_slice2 = det3(Gamma_slice2);
Gamma_slice2(1,1,:,:) = Gamma_slice2(1,1,:,:) + 1e-3;
cvals2 = getColorSPD(detGamma_slice2,N,cmap);
Gamma_slice2 = reshape(Gamma_slice2,3,3,N,1,N);
cvals2 = reshape(cvals2,N,1,N,3);

Gamma_slice3 = Gamma2(:,:,(Z == xx(floor(N/2))));
detGamma_slice3 = det3(Gamma_slice3);
Gamma_slice3(1,1,:,:) = Gamma_slice3(1,1,:,:) + 1e-3;
cvals3 = getColorSPD(detGamma_slice3,N,cmap);
Gamma_slice3 = reshape(Gamma_slice3,3,3,N,N);
% cvals3 = reshape(cvals3,N,N,3);

%%

if 1
gd = 1.;
f3 = figure(3);
plotSPD(gd*Gamma_slice1/max(Gamma_slice1(:)),'GridDistance',gd*[0,1,1],...
        'Colors',cvals1,'EllipsoidPoints',npts);
axis off;
view([90,0])
light('Position',[1 0 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;

print(f3,'aniso_tori_slice_1','-dpng','-r200');

%%
gd = 1.;
f3 = figure(3);
set(f3,'Position',[0,0,img_height,img_width]);
plotSPD(gd*Gamma_slice2/max(Gamma_slice2(:)),'GridDistance',gd*[1,0,1],...
    'Colors',cvals2,'EllipsoidPoints',npts);
axis off;
view([0,0])
light('Position',[0 -1 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;

print('aniso_tori_slice_2','-dpng','-r200');


%%
gd = 1.;
f3 = figure(3);
set(f3,'Position',[0,0,img_height,img_width]);
plotSPD(gd*Gamma_slice3/max(Gamma_slice3(:)),'GridDistance',gd*[1,1,0],...
    'Colors',cvals3,'EllipsoidPoints',npts);
axis off;
light('Position',[0 -1 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;
print('aniso_tori_slice_3','-dpng','-r200');

end
%%


N = 128;
g2 = myC('anisotori',20);
xx = linspace(-1,1,N);
[X,Y,Z] = ndgrid(xx,xx,xx);
Gamma2 = g2(X(:)',Y(:)',Z(:)');
Gamma2 = reshape(Gamma2,3,3,N^3);


tau = det3(Gamma2);
tau = reshape(tau,N,N,N);

f1 = figure(1); clf;
ezplot3D(tau,1);
colormap(cmap)
alpha(0.7)
adjustFigures(f1);
print(gcf,'tau_slices2.png','-dpng','-r150');


%%
N = 40;
g2 = myC('anisotori',20);
xx = linspace(-1,1,N);
[X,Y,Z] = ndgrid(xx,xx,xx);
Gamma2 = g2(X(:)',Y(:)',Z(:)');
Gamma2 = reshape(Gamma2,3,3,N^3);

Gamma_slice1 = Gamma2(:,:,(X == xx(floor(N/2))));
detGamma_slice1 = det3(Gamma_slice1);
cvals1 = getColorSPD(detGamma_slice1,N,cmap);
Gamma_slice1 = reshape(Gamma_slice1,3,3,1,N,N);
cvals1 = reshape(cvals1,1,N,N,3);

Gamma_slice2 = Gamma2(:,:,(Y == xx(floor(N/2))));
detGamma_slice2 = det3(Gamma_slice2);
Gamma_slice2(1,1,:,:) = Gamma_slice2(1,1,:,:) + 1e-3;
cvals2 = getColorSPD(detGamma_slice2,N,cmap);
Gamma_slice2 = reshape(Gamma_slice2,3,3,N,1,N);
cvals2 = reshape(cvals2,N,1,N,3);

Gamma_slice3 = Gamma2(:,:,(Z == xx(floor(N/2))));
detGamma_slice3 = det3(Gamma_slice3);
Gamma_slice3(1,1,:,:) = Gamma_slice3(1,1,:,:) + 1e-3;
cvals3 = getColorSPD(detGamma_slice3,N,cmap);
Gamma_slice3 = reshape(Gamma_slice3,3,3,N,N);
cvals3 = reshape(cvals3,N,N,3);

%%

if 1
gd = 1.;
f3 = figure(3);
set(f3,'Position',[0,0,img_height,img_width]);
plotSPD(gd*Gamma_slice1/max(Gamma_slice1(:)),'GridDistance',gd*[0,1,1],...
        'Colors',cvals1,'EllipsoidPoints',npts);
axis off;
view([90,0])
light('Position',[1 0 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;
print(f3,'aniso_tori2_slice_1','-dpng','-r200');

%%
gd = 1.;
f3 = figure(3);
set(f3,'Position',[0,0,img_height,img_width]);
plotSPD(gd*Gamma_slice2/max(Gamma_slice2(:)),'GridDistance',gd*[1,0,1],...
    'Colors',cvals2,'EllipsoidPoints',npts);
axis off;
view([0,0])
light('Position',[0 -1 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;
print('aniso_tori2_slice_2','-dpng','-r200');


%%
gd = 1;
f3 = figure(3);
set(f3,'Position',[0,0,img_height,img_width]);
plotSPD(gd*Gamma_slice3/max(Gamma_slice3(:)),'GridDistance',gd*[1,1,0],...
    'Colors',cvals3,'EllipsoidPoints',npts);
axis off;
light('Position',[0 -1 0],'Style','infinite');

set(f3,'Position',[0,0,img_width,img_height]);
f3.Children(1).Position(1) = 0.05;
f3.Children(1).Position(2) = 0.05;
f3.Children(1).Position(3) = 0.9;
f3.Children(1).Position(4) = 0.9;
print('aniso_tori2_slice_3','-dpng','-r200');
end
