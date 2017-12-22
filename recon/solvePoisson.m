function [x,A] = solvePoisson(x,N,varargin)
% x is 3d object (of size N^3 x 1) with proper boundary conditions


x = x(:);           % flatten x
h = 2/(N-1);        % grid-width
tol = 1e-5;         % tolerance
maxit = 3000;       % max number of iterations

disp('building stiffness matrix...'); 
tic
D = spdiags(ones(N,1)/(2*h),1,N,N) + spdiags(-ones(N,1)/(2*h),-1,N,N);
D(1,1:3) =     [-3  4 -1]/(2*h);
D(N,(N-2):N) = [ 1 -4  3]/(2*h);

I = speye(N);
Dx = kron(I,kron(I,D));
Dy = kron(I,kron(D,I));
Dz = kron(D,kron(I,I));

A = (Dx*Dx + Dy*Dy + Dz*Dz);

disp('Stiffness matrix done.')
toc;

if ~isempty(varargin) 
   b = varargin{1};
   b = b(:);
else
   b = zeros(N^3,1);
end

b = b - A*x;

% get logical indices for boundaries
[xbd0,xbd1,ybd0,ybd1,zbd0,zbd1] = getBdry(N);
J = ones(N^3,1);
J([xbd0,xbd1,ybd0,ybd1,zbd0,zbd1]) = 0;
J = logical(J);

disp('Solving system...')
tic; 
x(J) = bicgstab(A(J,J),b(J),tol,maxit); 
toc;

end
