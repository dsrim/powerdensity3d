function [dAdx,dAdy,dAdz] = computeGradMat(A,h)
%COMPUTEGRADMAT
%
% compute gradients of a matrix function over 3D domain

M = size(A,3);          % # of 3D gridpts
N = floor(M^(1/3));     % # of grid pts in each 1D axis
if (M ~= N^3)           % correct for possible error
    N = N+1;
end


dAdx = zeros(3,3,M);
dAdy = zeros(3,3,M);
dAdz = zeros(3,3,M);

for i = 1:3
    for j = 1:3

        B = squeeze(A(i,j,:));
        
        [dBdx, dBdy, dBdz] = computeGrad(B,h);
        
        dAdx(i,j,:) = dBdx;
        dAdy(i,j,:) = dBdy;
        dAdz(i,j,:) = dBdz;
    end
end


end