function C = outer3(A,B)
%OUTER3
%
% compute outer product with extra index

M = size(A,2);

C = zeros(3,3,M);

for i = 1:3
  for j = 1:3
    C(i,j,:) = A(i,:).*B(j,:);
  end
end

end