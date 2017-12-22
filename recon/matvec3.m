function C = matvec3(A,B)
%MATMUL3
%
% do matrix multiplcation for 3x3 matrix and 3-vector, with extra index
% the extra index should be huge

C = zeros(size(B));

for i = 1:3
    C(i,:) = squeeze(A(i,1,:))'.*B(1,:) ...
           + squeeze(A(i,2,:))'.*B(2,:) ...
           + squeeze(A(i,3,:))'.*B(3,:);
end


end