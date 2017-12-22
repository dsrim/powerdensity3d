function G = computeGStab(B,H)

M = size(H,3);
G = zeros(3,3,M);

G(:,1,:) = matvec3(B, cofactor3(H, squeeze(B(1,:,:))));
G(:,2,:) = matvec3(B, cofactor3(H, squeeze(B(2,:,:))));
G(:,3,:) = matvec3(B, cofactor3(H, squeeze(B(3,:,:))));

end