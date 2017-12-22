function A = convertAij(A11,A12,A13,A21,A22,A23,A31,A32,A33)
%CONVERTAIJ
%
% convert matrix format

M = numel(A11);
A = zeros(3,3,M);

A(1,1,:) = A11(:);
A(1,2,:) = A12(:);
A(1,3,:) = A13(:);

A(2,1,:) = A21(:);
A(2,2,:) = A22(:);
A(2,3,:) = A23(:);

A(3,1,:) = A31(:);
A(3,2,:) = A32(:);
A(3,3,:) = A33(:);

end
