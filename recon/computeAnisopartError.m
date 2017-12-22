function [e4tgamma, Tgamma] = computeAnisopartError(Gamma,tgamma)

M = size(Gamma,3);
Tau = det3(Gamma).^(1/3);
Tgamma = Gamma./reshape(repmat(Tau,9,1),3,3,M);       % true

detTgamma = det3(tgamma).^(1/3);
tgamma = tgamma./reshape(repmat(detTgamma,9,1),3,3,M);   % approx


e4tgamma = Tgamma - tgamma;

end
