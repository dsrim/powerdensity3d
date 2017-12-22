function detA = det3(A)
%DET3
%
% Compute det(A) for 3x3 matrix A with long extra index

A11 = squeeze(A(1,1,:))';
A12 = squeeze(A(1,2,:))';
A13 = squeeze(A(1,3,:))';

A21 = squeeze(A(2,1,:))';
A22 = squeeze(A(2,2,:))';
A23 = squeeze(A(2,3,:))';

A31 = squeeze(A(3,1,:))';
A32 = squeeze(A(3,2,:))';
A33 = squeeze(A(3,3,:))';

detA =   A11.*(A22.*A33 - A23.*A32) ...
       - A12.*(A21.*A33 - A23.*A31) ...
       + A13.*(A21.*A32 - A22.*A31);

end
