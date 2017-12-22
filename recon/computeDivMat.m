function divA = computeDivMat(A,h)

N = size(A.m11,1);
M = numel(A.m11);

divA = initMatx(3,3,N);


for i = 1:3
    for j = 1:3
        istr = num2str(i);
        jstr = num2str(j);
        B = getMatrixIJ(A,i,j);
        
        [dBdx, dBdy, dBdz] = computeGrad(B,h);
        
        eval(['divA.m' istr jstr ' = dBdx + dBdy + dBdz;'])
    end
end



end