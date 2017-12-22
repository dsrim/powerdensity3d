function Mu = computeAnisoMu(H)
%COMPUTEANSIOMU
%
% compute coordinates \mu by solving H(1:3,1:3) \ H(4:,:)
%

K = size(H,1) - 3;
M = size(H,3);
Mu = zeros(3,K,M);

for k = 1:K
   Mu(:,k,:) = cramersrule3(H(1:3,1:3,:),squeeze(H(1:3,k+3,:)));
end


end
