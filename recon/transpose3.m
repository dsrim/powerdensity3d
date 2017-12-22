function B = transpose3(A)
%TRANSPOSE3
%
% do matrix transpose for 3x3 matrices, with extra index
% the extra index should be huge

n1 = size(A,1);
n2 = size(A,2);
M = size(A,3);
B = zeros(n1,n2,M);

for i = 1:n1
    for j = 1:n2
        B(j,i,:) = A(i,j,:);
    end
end

end
