function [u,sigval] = solveScalarConductivity(u,b,sigma_type,const)
%SOLVESCALARCONDUCTIVITY
%
% solve the scalar conductivity problem for given sigma on a cube:
%
%     \div \sigma \grad u = b
%                       u = g
%
% using simple finite difference
%
% u:    solution array of dim N x N x N, zero for interior pts
%       non-zero entries for Dirichlet boundary info.
% b:    source array of dim N x N x N, 
% sigmatype:    type of sigma to choose (help iso_conductivity)
% const:        parameter for sigmatype (help iso_conductivity)
%
% NOTE: the output solution is in meshgrid ordering

w0 = -1;            % left boundary 
w1 =  1;            % right boundary
N = size(u,1);      % get N from argument u
h = 2/(N-1);        % grid width
maxit = 2000;       % max iteration for biCG
tol = 1e-8;         % tolerance for biCG

x = linspace(w0,w1,N);  % 1D grid

[ X, Y, Z] = meshgrid(x,x,x); 
[iX,iY,iZ] = meshgrid(1:N,1:N,1:N);

iX1 = iX(:); 
iY1 = iY(:); 
iZ1 = iZ(:);

% evaluate sigma at grid-points
[sigval,sigfun] = iso_conductivity(sigma_type,const,[X(:),Y(:),Z(:)]');

% get logical indices for boundary
n4Ds = getBdryAll(N);

tic;
disp(['(' mfilename ') building matrix']);

Ns = 1:N^3;
M = N^3;
Aloc = zeros(5,3,M);
id1  =  zeros(15,M);
id2  =  zeros(15,M);

for j = Ns
  
  if n4Ds(j,1)
      % lies on the bdry, nothing to do.
  else
    sr = n4Ds(j,2:4)';
    sl = n4Ds(j,5:7)';
    w = [sr, sl];
    ii = (repmat([N;1;N^2],1,5).*(repmat(-2:2,3,1) ...     % centered stencil
                 + repmat(sr,1,5) - repmat(sl,1,5)))';     % shift stencil 
    I =  ii(:) + j;
    id2(:,j) = I;
    cs = [X(I),Y(I),Z(I)];
    xp = .5*(cs([1:4,6:9,11:14],:) + cs([2:5,7:10,12:15],:));
    sig = reshape(sigfun(xp'),4,3);
    
    for k = 1:3
      if w(k,1) 
        Aloc(:,k,j) = [11/12,    0,    0,     0;
                      -11/12, -3/4,    0,     0;
                           0,  3/4, -1/4,     0;
                           0,    0,  1/4,  1/12;
                           0,    0,    0, -1/12]/h^2*sig(:,k);
      elseif w(k,2)
        Aloc(:,k,j) =  [-1/12,    0,    0,    0;
                         1/12,  1/4,    0,    0; 
                            0, -1/4,  3/4,    0;
                            0,    0, -3/4,-11/12;
                            0,    0,    0, 11/12]/h^2*sig(:,k); 
      else
        Aloc(:,k,j) = [-1/12,   0,     0,     0; 
                        1/12, 5/4,     0,     0; 
                           0,-5/4,  -5/4,     0;
                           0,   0,   5/4,  1/12; 
                           0,   0,     0, -1/12]/h^2*sig(:,k);
      end
      
    end
    id1(:,j) = j*ones(15,1);
          
  end
end

J = n4Ds(:,1);
Aloc = Aloc(:,:,~J);
A = sparse(id1(:,~J),id2(:,~J),Aloc(:),M,M); % write stiffness matrix
toc;

b = b(:);           % reshape RHS
b = b - A*u(:);     % apply BC

disp(['(' mfilename ') Solving one conductivity problem...'])
tic; 
u(~J) = bicg(A(~J,~J),b(~J),tol,maxit); 
toc;

u = reshape(u,N,N,N);       % reshape solution

end