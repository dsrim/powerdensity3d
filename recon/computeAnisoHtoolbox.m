function H = computeAnisoHtoolbox(Gamma,h,N,varargin)
%COMPUTEANISOHTOOLBOX
%
% compute anisotropic power density matrix H for arbitrary number of 
% solutions
%

M = N^3;
K = length(varargin);

xx = linspace(-1,1,N);
[X,Y,Z] = meshgrid(xx,xx,xx);

% note: ndgrid ordering
for k = 1:K
    eval(['u' num2str(k) ' = varargin{' num2str(k) '};']);
    eval(['[du' num2str(k) 'dx,du' num2str(k) 'dy,du' num2str(k) 'dz] ' ...
          '= u' num2str(k) '.evaluateGradient(X,Y,Z);']);
    eval(['du' num2str(k) 'dx = reshape(du' num2str(k) 'dx,N,N,N);'])
    eval(['du' num2str(k) 'dy = reshape(du' num2str(k) 'dy,N,N,N);'])
    eval(['du' num2str(k) 'dz = reshape(du' num2str(k) 'dz,N,N,N);'])

    eval(['du' num2str(k) 'dx = permute(du' num2str(k) 'dx,[2 1 3]);'])
    eval(['du' num2str(k) 'dy = permute(du' num2str(k) 'dy,[2 1 3]);'])
    eval(['du' num2str(k) 'dz = permute(du' num2str(k) 'dz,[2 1 3]);'])
end


H = zeros(K,K,M);

for i = 1:K
  for j = 1:K

    eval(['H(' num2str(i) ',' num2str(j) ',:) = ' ...
 '( squeeze(Gamma(1,1,:)).* du' num2str(i) 'dx(:).*du' num2str(j) 'dx(:)' ...
 '+ squeeze(Gamma(1,2,:)).* du' num2str(i) 'dx(:).*du' num2str(j) 'dy(:)' ...
 '+ squeeze(Gamma(1,3,:)).* du' num2str(i) 'dx(:).*du' num2str(j) 'dz(:)' ...
 '+ squeeze(Gamma(2,1,:)).* du' num2str(i) 'dy(:).*du' num2str(j) 'dx(:)' ...
 '+ squeeze(Gamma(2,2,:)).* du' num2str(i) 'dy(:).*du' num2str(j) 'dy(:)' ...
 '+ squeeze(Gamma(2,3,:)).* du' num2str(i) 'dy(:).*du' num2str(j) 'dz(:)' ...
 '+ squeeze(Gamma(3,1,:)).* du' num2str(i) 'dz(:).*du' num2str(j) 'dx(:)' ...
 '+ squeeze(Gamma(3,2,:)).* du' num2str(i) 'dz(:).*du' num2str(j) 'dy(:)' ...
 '+ squeeze(Gamma(3,3,:)).* du' num2str(i) 'dz(:).*du' num2str(j) 'dz(:) );'])

  end
end

end
