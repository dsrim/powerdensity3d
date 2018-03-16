function [B,detvalH,detvalB,condA,wflag] = computePerpB(Z,ZHOm,H)
%COMPUTEPERPB
%
% returns array B of size 3 x 3 x M
%       where M = N^3 is the number of 3D grid points
%
% each 3x3 matrix B(:,:,j) in vectorized form, 
% i.e. vector(B(:,:,j))), is orthogonal to the matrices
%       vector(Z(:,:,k,j)), vector(ZHOm(:,:,l,k,j))
%
% the orthogonal matrix is found by SVD

M = size(H,3);          % # of 3D gridpts
N = floor(M^(1/3));     % # of grid pts in each 1D axis
if (M ~= N^3)           % correct for possible error
    N = N+1;
end

wflag = zeros(M,1);              % flag warning

K = size(Z,3);

B = zeros(3,3,M);
detvalH = zeros(M,1);
detvalB = zeros(M,1);
condA = zeros(M,1);

display(['(computePerpB) number of redundant sols: ' num2str(K)])

parfor j = 1:M
    A = zeros(9,K + 3*K);   % matrix to contain vec(Z)s and vec(ZHOms)s as cols

    Hx = H(:,:,j);
    Hx = Hx(1:3,1:3);
    
    for k = 1:K
        A(:,k) = reshape(Z(:,:,k,j),9,1);
        ZHOmx = ZHOm(:,:,:,k,j);
        A(:,(3*k):(3*k + 2)) = reshape(ZHOmx(:,:,1:3),9,3);
    end
    
    %[Q,~] = qr(A);
    [U,S,~] = svd(A);

    tol = max(size(A(:,:))) * eps(norm(A(:,:)));
    if (sum(diag(S) > tol) < 8)     % if rank of A is less than 8
%         display(['warning: rank(A) = ' num2str(rank(A))])
        wflag(j) = 1;
    end
	
 	Bx = U(:,9);
	Bx = U(:,9)/det(reshape(Bx,3,3));
	for l = 1:10
	for k = 1:size(A,2)
		w = A(:,k)/norm(A(:,k));
	    Bx = Bx - w*w'*Bx;
	    Bx = Bx/abs(det(reshape(Bx,3,3)))^(1/3);
	end
	end
    Bx = reshape(Bx,3,3);

    detBx = computeCubicsqrt(sqrt(det(Hx))/det(Bx));
    
    if abs(detBx) > 1e1
        %display('detBx > 1e1')
        wflag(j) = 2;
    end
    Bx = detBx*Bx;
    
    B(:,:,j) = Bx;
    detvalH(j) = det(Hx);
    detvalB(j) = detBx;
    condA(j) = S(1,1)/S(K,K);

end

end




function val = computeCubicsqrt(x)

    if x >= 0
        val = x^(1/3);
    else
        val = - abs(x)^(1/3);
    end

end



