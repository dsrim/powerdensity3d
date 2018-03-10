function computeIsoSols(sigma_type,const,N,hmax)
%COMPUTEISOSOLS
%
% Generate three solutions and related functions (gradients, S, detS, etc.)
% the isotropic (scalar) conductivity problems. 
% Each solution corresponding to simple boundary conditions x,y,z.
%
% saves output filename IsoSols_[sigma_type]_[const]_[N].mat

% save soutions to save_fname
save_fname = getSaveFname('IsoSols',sigma_type,const,N);

%% generate 3 solutions


% old FDM code
if 0
h = 2 / (N-1);
w0 = -1; 
w1 =  1;
x = linspace(w0,w1,N);
[X,Y,Z] = meshgrid(x,x,x);

% initialize solution variables
u1 = zeros(N,N,N);
u2 = zeros(N,N,N);
u3 = zeros(N,N,N);

% set boundary conditions
J = getBdryAll(N);
J = J(:,1);
u1(J) = X(J);
u2(J) = Y(J);
u3(J) = Z(J);

b = zeros(N,N,N);    % RHS is zero

% solve three conductivity problems
[u1,sigval] = solveScalarConductivity(u1,b,sigma_type,const);
[u2,     ~] = solveScalarConductivity(u2,b,sigma_type,const);
[u3,     ~] = solveScalarConductivity(u3,b,sigma_type,const);

sigval = reshape(sigval,N,N,N);
end

% get function handle for isotropic conductivity 
[~,cscalar_fctn] = iso_conductivity(sigma_type, const,[1,1,1]'); 
c = @(r,s) [1;0;0;0;1;0;0;0;1]*cscalar_fctn([r.x(:),r.y(:),r.z(:)]');  

%% Set up PDE Toolbox solver
%
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

%% interpolate FE solution to uniform grid
x = linspace(-1,1,N);
h = 2/(N-1);
[X,Y,Z] = meshgrid(x,x,x);
sigval  = reshape(cscalar_fctn([X(:),Y(:),Z(:)]'),N,N,N); 

u1 = result1.interpolateSolution([X(:),Y(:),Z(:)]');
u2 = result2.interpolateSolution([X(:),Y(:),Z(:)]');
u3 = result3.interpolateSolution([X(:),Y(:),Z(:)]');

u1 = reshape(u1,N,N,N);
u2 = reshape(u2,N,N,N);
u3 = reshape(u3,N,N,N);


%% compute gradients 
disp(['(' mfilename ') computing gradients'])
[du1dx, du1dy, du1dz] = computeGrad(u1,h,'meshgrid');
[du2dx, du2dy, du2dz] = computeGrad(u2,h,'meshgrid');
[du3dx, du3dy, du3dz] = computeGrad(u3,h,'meshgrid');
disp(['(' mfilename ') done.'])


%% compute detS, detDU
disp(['(' mfilename ') computing detS, detDU'])
detS = zeros(N^3,1);
detDU = zeros(N^3,1);

for k = 1:N^3
  S = sqrt(sigval(k)).*[du1dx(k), du1dy(k), du1dz(k);
                        du2dx(k), du2dy(k), du2dz(k);
                        du3dx(k), du3dy(k), du3dz(k)]';
  DU = [du1dx(k), du1dy(k), du1dz(k);
        du2dx(k), du2dy(k), du2dz(k);
        du3dx(k), du3dy(k), du3dz(k)]';
  detS(k,1)  = det(S);
  detDU(k,1) = det(DU);
end
disp(['(' mfilename ') done.'])

[mindet, ~] = min(detS);
[maxdet, ~] = max(detS);
display(['min/max det(S): ' num2str(mindet) ' / ' num2str(maxdet)])

%% Save outputs
disp(['(' mfilename ') plotting / saving outputs'])
save(save_fname,'u1','u2','u3','detS','X','Y','Z','h','detDU','sigval',...
                'du1dx','du2dx','du3dx',...
                'du1dy','du2dy','du3dy',...
                'du1dz','du2dz','du3dz','-v7.3')

end

