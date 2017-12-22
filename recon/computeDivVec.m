function divA = computeDivVec(A,h)

if size(A,1) > size(A,2)
  A = A';
end

B1 = squeeze(A(1,:));
B2 = squeeze(A(2,:));
B3 = squeeze(A(3,:));
        
[dB1dx, ~, ~] = computeGrad(B1,h);
[~, dB2dy, ~] = computeGrad(B2,h);
[~, ~, dB3dz] = computeGrad(B3,h);
        
divA = dB1dx + dB2dy + dB3dz;



end