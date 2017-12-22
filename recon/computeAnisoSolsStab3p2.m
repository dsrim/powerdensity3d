function computeAnisoSolsStab3p2(gamma_type,const,N,hmax)
%COMPUTEANISOSOLS3P2
%
% compute solutions to the anisotropic conductivity problem
% for the boundary conditions x,y,z,(x+2),(y+2) for the 3+2 algorithm
% using MATLAB PDE Toolbox.
%
% (tested on MATLAB R2016b)
             
% get function handle for anisotropic conductivity Gamma
cfunction = aniso_conductivity(gamma_type, const);
c = @(r,s) cfunction(r.x,r.y,r.z);

% set up uniform Cartesian grid of size N x N x N
M = N^3;
xx = linspace(-1,1,N);
h = 2/(N-1);                        % grid width
[X,Y,Z] = ndgrid(xx,xx,xx);         % 3D grid

% evaluate Gamma on uniform grid
r0.x = X(:)';  
r0.y = Y(:)';  
r0.z = Z(:)';


%% Set up PDE Toolbox solver

% Create 3-D copies of the remaining mesh grid points, with the
% |z|-coordinates ranging from 0 through 1. Together, these points
% constitute a 3-D point cloud. Combine the points into an |alphaShape|
% object.
[xg, yg] = meshgrid(-1:.25:1);
xg = xg(:);
yg = yg(:);
zg = ones(numel(xg),1);
xg = repmat(xg,9,1);
yg = repmat(yg,9,1);
zg = zg*(-1:.25:1);
zg = zg(:);
shp = alphaShape(xg,yg,zg);

% Obtain a surface mesh of the |alphaShape| object.
[elements,nodes] = boundaryFacets(shp);

% Put the data in the correct shape for |geometryFromMesh|.
nodes = nodes';
elements = elements';

% Create a PDE model and import the surface mesh.
model = createpde();
geometryFromMesh(model,nodes,elements);

% generate the mesh (last arg = max edge length)
generateMesh(model,'Hmax',hmax)


% set PDE coefficients
f = 0;      % no source term
a = 0;      % no advection term
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

%% solve PDEs
% bc x
bc = @(r,s) r.x;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result1 = solvepde(model);
toc;

% bc y
bc = @(r,s) r.y;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result2 = solvepde(model);
toc;

% bc z
bc = @(r,s) r.z;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result3 = solvepde(model);
toc;

% bc (x+2)(y+2)
bc = @(r,s) (r.x + 2).*(r.y + 2);
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result4 = solvepde(model);
toc;

% bc (y+2)(z+2)
bc = @(r,s) (r.y + 2).*(r.z + 2);
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result5 = solvepde(model);
toc;

% bc (z+2)(x+2)
bc = @(r,s) (r.z + 2).*(r.x + 2);
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result6 = solvepde(model); 
toc;

% bc x+2.5*(z+1)^2
bc = @(r,s) r.x + 2.5*(r.z + 2).^2;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result21 = solvepde(model); 
toc;

% bc y+2.5*(x+1)^2
bc = @(r,s) r.y + 2.5*(r.x + 2).^2;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result22 = solvepde(model); 
toc;

% bc z+2.5*(y+1)^2
bc = @(r,s) r.z + 2.5*(r.y + 2).^2;
applyBoundaryCondition(model,'dirichlet','Face',1:6,'u',bc);
disp('solving for bc:'); display(bc); 
tic;
result23 = solvepde(model); 
toc;

save_fname = getSaveFname('AnisoSolsStab3p2',gamma_type,const,N);
save(save_fname, ...
      'result1', 'result2', 'result3', 'result4', 'result5','result6',...
     'result21','result22','result23');

end