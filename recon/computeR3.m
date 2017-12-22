function R = computeR3(X)
%COMPUTER3

orig_size = size(X);
M = size(X(1,:),2);      % number of grid points on the contour integral
R = zeros(3,3,M);

%  R = [aa^2 + bb^2 - cc^2 - dd^2, 2*(bb*cc-aa*dd), 2*(bb*dd+aa*cc); ...
%      2*(bb*cc + aa*dd), aa^2 - bb^2 + cc^2 - dd^2, 2*(cc*dd-aa*bb); ...
%      2*(bb*dd-aa*cc), 2*(cc*dd+aa*bb), aa^2 - bb^2 - cc^2 + dd^2];

R(1,1,:) = X(1,:).^2 + X(2,:).^2 - X(3,:).^2 - X(4,:).^2;
R(2,2,:) = X(1,:).^2 - X(2,:).^2 + X(3,:).^2 - X(4,:).^2;
R(3,3,:) = X(1,:).^2 - X(2,:).^2 - X(3,:).^2 + X(4,:).^2;

R(1,2,:) = 2*(X(2,:).*X(3,:) - X(1,:).*X(4,:));
R(2,1,:) = 2*(X(2,:).*X(3,:) + X(1,:).*X(4,:));

R(1,3,:) = 2*(X(2,:).*X(4,:) + X(1,:).*X(3,:));
R(3,1,:) = 2*(X(2,:).*X(4,:) - X(1,:).*X(3,:));

R(2,3,:) = 2*(X(3,:).*X(4,:) - X(1,:).*X(2,:));
R(3,2,:) = 2*(X(3,:).*X(4,:) + X(1,:).*X(2,:));

% reshape to get original size
R = reshape(R,[3,3,orig_size(2:end)]);

end
