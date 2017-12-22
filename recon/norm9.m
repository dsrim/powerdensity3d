function normA = norm9(A)
%NORM3
%
% compute l2 norm of 3x3xM array by flattening the 3x3 matrix as
% a 9-vector

M = size(A,3);
A = reshape(A,9,M);

normA = zeros(M,1);
for j = 1:M
    normA(j) = norm(A(:,j));
end

end