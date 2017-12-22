function C = matmul3(A,B)
%MATMUL3
%
% do matrix multiplcation for 3x3 matrices, with extra index
% the extra index should be huge

C = zeros(size(A));

for i = 1:3
  for j = 1:3
    C(i,j,:) = A(i,1,:).*B(1,j,:) + A(i,2,:).*B(2,j,:) + A(i,3,:).*B(3,j,:);
  end
end


end