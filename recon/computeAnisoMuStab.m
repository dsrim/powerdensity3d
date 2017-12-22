function tMu = computeAnisoMuStab(H)
%COMPUTEANSIOMUSTAB
%
% compute coordinates \mu by solving H(1:3,1:3) \ H(4:,:)
% 

K = size(H,1) - 3;
M = size(H,3);
tMu = zeros(3,K,M);

for k = 1:K
   tMu(:,k,:) = cofactor3(H(1:3,1:3,:), squeeze(H(1:3,k+3,:)));
end


end
