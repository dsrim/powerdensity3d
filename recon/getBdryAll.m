function n4Ds = getBdryAll(N)
%GETBDRYALL
%
% get logical indices for the boundary considering it as 
% a N x N x N meshgrid array flattened to be of dim N^3 x 1
%
% N:    grid-size along one axis (total grid-size equals N^3)

[iX,iY,iZ] = meshgrid(1:N,1:N,1:N);
iX1 = iX(:); iY1 = iY(:); iZ1 = iZ(:);

n4Ds = false(N^3,7);

n4Ds(((iX1 == 1) + (iX1 == N) ...
    + (iY1 == 1) + (iY1 == N) ...
    + (iZ1 == 1) + (iZ1 == N)) > 0,1) = true;
n4Ds((iX1 == 2),  2) = true;
n4Ds((iY1 == 2),  3) = true;
n4Ds((iZ1 == 2),  4) = true;
n4Ds((iX1 == N-1),5) = true;
n4Ds((iY1 == N-1),6) = true;
n4Ds((iZ1 == N-1),7) = true;

end