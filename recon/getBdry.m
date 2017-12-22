function [xbd0,xbd1,ybd0,ybd1,zbd0,zbd1] = getBdry(N)
%GETBDRY
%
% get boundary indices for a cubic domain
%
% N:    grid-size along one axis (total grid-size equals N^3)
    
xbd0 = 1:N^2; xbd1 = (N^3 - N^2+1):N^3;
ybd0 = 1:N:N^3; ybd1 = N:N:N^3;
zbd0 = [];
zbd1 = [];
for j = 1:N
    zbd0 = [zbd0, j:N^2:N^3]; 
    zbd1 = [zbd1, (N^2 - N + j):N^2:N^3];
end
end