function [e4tgamma,e4gamma,e4tau] = computeAnisoError(Gamma,tgamma,tau)
%COMPUTEANISOERROR

M = numel(tau);
Tau = det3(Gamma).^(1/3);

[e4tgamma,Tgamma] = computeAnisopartError(Gamma,tgamma);

e4tau = Tau(:) - tau(:);
e4gamma = Gamma - tgamma.*reshape(repmat(tau(:)',9,1),3,3,M);

end
