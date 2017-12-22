function cvals = getColorSPD(V,N,cmap)

Nc = size(cmap,1);

minV = min(V(:));
maxV = max(V(:));
V = (V - minV)*(1 - 1e-2) + minV + 1e-2;
cX = linspace(minV,maxV,Nc);
cvals = zeros(N,N,3);
for j = 1:3
    cvals(:,:,j) = reshape(interp1(cX',cmap(:,j),V,'linear'),N,N);
end

end