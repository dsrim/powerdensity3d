function C = cofactor3(A,b)
%COFACTOR3
%
% Compute cofactors for 3x3 system with long extra index

M = size(A,3);
C = zeros(3,M);

A11 = squeeze(A(1,1,:))';
A12 = squeeze(A(1,2,:))';
A13 = squeeze(A(1,3,:))';

A21 = squeeze(A(2,1,:))';
A22 = squeeze(A(2,2,:))';
A23 = squeeze(A(2,3,:))';

A31 = squeeze(A(3,1,:))';
A32 = squeeze(A(3,2,:))';
A33 = squeeze(A(3,3,:))';

b1 = squeeze(b(1,:));
b2 = squeeze(b(2,:));
b3 = squeeze(b(3,:));


C(1,:) = ( b1.*(A22.*A33 - A23.*A32) ...
        - A12.*(b2 .*A33 - A23.*b3 ) ...
        + A13.*(b2 .*A32 - A22.*b3 ));

C(2,:) = (A11.*(b2 .*A33 - A23.*b3 ) ...
        - b1 .*(A21.*A33 - A23.*A31) ...
        + A13.*(A21.*b3  - b2 .*A31));
    
C(3,:) = (A11.*(A22.*b3  - b2 .*A32) ...
        - A12.*(A21.*b3  - b2 .*A31) ...
        + b1 .*(A21.*A32 - A22.*A31));

end
