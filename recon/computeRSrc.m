function b = computeRSrc(tgamma,H,B)

M = size(H,3);          % # of 3D gridpts
N = floor(M^(1/3));     % # of grid pts in each 1D axis
if (M ~= N^3)           % correct for possible error
    N = N+1;
end
h = 2/(N-1); % bad

b = zeros(3,M);         % RHS vector
iH = zeros(3,3,M);      % inverse of H
I = getGamma(N,'Id');   % the identity

for j = 1:3
  iH(:,j,:) = cramersrule3(H,squeeze(I(:,j,:))); 
end

ldH = log(det3(H));     % log(det(H))

[dldHdx,dldHdy,dldHdz] = computeGrad(ldH,h);
[diHdx,diHdy,diHdz] = computeGradMat(iH,h);

C =  matvec3(B, matvec3(diHdx,squeeze(B(1,:,:))) ...
              + matvec3(diHdy,squeeze(B(2,:,:))) ...
              + matvec3(diHdz,squeeze(B(3,:,:))));

b = 1/3 * matvec3(tgamma, [dldHdx'; dldHdy'; dldHdz']) + 2/3 * C;       
b = b';   % kluging it

end
