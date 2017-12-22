function normA = norm3(A)
%NORM3
%
% compute l2 norm with extra index

normA = sqrt(A(1,:).^2 + A(2,:).^2 + A(3,:).^2);

end