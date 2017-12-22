function computeIsoSols(sigma_type,const,N)
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
                'du1dz','du2dz','du3dz')

end

