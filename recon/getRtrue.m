function Rtrue = getRtrue(sigma_type,const,N)
%GETRTRUE
%
% helper function that loads the "true" R matrix
% requires output from computeTij

Tij_fname = getSaveFname('Tij',sigma_type,const,N);
load(Tij_fname)

Rtrue = zeros(3,3,N^3);

Rtrue(1,1,:) = R11(:);
Rtrue(1,2,:) = R12(:);
Rtrue(1,3,:) = R13(:);

Rtrue(2,1,:) = R21(:);
Rtrue(2,2,:) = R22(:);
Rtrue(2,3,:) = R23(:);

Rtrue(3,1,:) = R31(:);
Rtrue(3,2,:) = R32(:);
Rtrue(3,3,:) = R33(:);

end
