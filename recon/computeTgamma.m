function tgamma = computeTgamma(B,H)

M = size(H,3);
N = M^(1/3);
if (M ~= N^3)
    N = N + 1;
end

tgamma = zeros(3,3,M);


tgamma(:,1,:) = matvec3(B, cramersrule3(H, squeeze(B(1,:,:))));
tgamma(:,2,:) = matvec3(B, cramersrule3(H, squeeze(B(2,:,:))));
tgamma(:,3,:) = matvec3(B, cramersrule3(H, squeeze(B(3,:,:))));


end
