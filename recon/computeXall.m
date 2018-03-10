function computeXall(sigma_type,const,N)
%COMPUTEXALL
%
% compute quaternionic dynamics system along all curves

X_fname = getSaveFname('X',sigma_type,const,N);

Xall = zeros(4,N,N,N);          % initialize variable

H = getH(sigma_type,const,N);         % obtain power density data
Vijl = getVijl(sigma_type,const,N);   % compute/load Vijl
Ball = computeB(H);                   % compute matrix B (no dep. on q)

parfor j = 1:N
  for k = 1:N
      
    % set curve
    curve = zeros(N,4);

    % single index N*(i-1) + j + N^2*(k-1) for grid point (i,j,k)
    % (meshgrid index ordering)
    % integrate along x-direction 
    curve(:,1) = N*((1:N) - 1) + j + N^2*(k - 1);    

    % x-gradient                                  
    curve(:,2) = 1; 

    % solve system for one curve
    Xall(:,j,:,k) = computeX(curve,Vijl,Ball);

  end
end
fprintf('\n = done.')

% save output
save(X_fname, 'Xall','-v7.3')

end
