function b = computeRSrcStab(G,H,B)

M = size(H,3);          % # of 3D gridpts
N = floor(M^(1/3));     % # of grid pts in each 1D axis
if (M ~= N^3)           % correct for possible error
    N = N+1;
end
h = 2/(N-1); % bad

b = zeros(3,M);         % RHS vector
I = getGamma(N,'Id');   % the identity
detH = det3(H);         % |H|

for j = 1:3
  tH(:,j,:) = cofactor3(H(1:3,1:3,:),squeeze(I(:,j,:))); 
end

[ddetHdx,ddetHdy,ddetHdz] = computeGrad(detH,h);
[dtHdx,dtHdy,dtHdz] = computeGradMat(tH,h);

C =  repmat(detH,3,1).*matvec3(B, matvec3(dtHdx,squeeze(B(1,:,:))) ...
                                + matvec3(dtHdy,squeeze(B(2,:,:))) ...
                                + matvec3(dtHdz,squeeze(B(3,:,:))));

b = - 1/3 * matvec3(G, [ddetHdx'; ddetHdy'; ddetHdz']) + 2/3 * C;       
b = b';   % kluging it

end
