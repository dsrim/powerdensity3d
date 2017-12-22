function Mu = computeAnisoMuPermute(H, sInd)
%COMPUTEANSIOMUPERMUTE
%
% compute coordinates \mu by solving H(1:3,1:3) \ H(4:,:)
%

K = size(H,1) - 3;
M = size(H,3);

Mu = zeros(3,K,M);

parfor j = 1:M
  sIndx = sInd(:,j);
  Hx = H(sIndx(1:3), sIndx(1:3), j);
  Hb = H(sIndx(1:3), sIndx(4:5), j);
  Mu(:,:,j) = Hx \ Hb;
end


end
