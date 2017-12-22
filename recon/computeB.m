function B = computeB(H)
%COMPUTEB
%
% computes B that depends only on log(det(H))
% B needs only be computed once for all spatial points

N = size(H,3);          % number of grid pts (dof)
M = N^3;
h = 2/(N-1);            % bad

H = reshape(H,3,3,N^3);
ldH = log(det3(H));     % log(det(H))
ldH = reshape(ldH,N,N,N);

[dldHdx,dldHdy,dldHdz] = computeGrad(ldH,h,'meshgrid');

B = zeros(4,4,M,3);

% m = 1
B(1,3,:,1) = - dldHdz(:);
B(2,4,:,1) =   dldHdz(:);

B(1,4,:,1) =   dldHdy(:);
B(2,3,:,1) =   dldHdy(:);

B(3,1,:,1) = - B(1,3,:,1);
B(4,2,:,1) = - B(2,4,:,1);

B(4,1,:,1) = - B(1,4,:,1);
B(3,2,:,1) = - B(2,3,:,1);


% m = 2
B(1,2,:,2) = + dldHdz(:);
B(3,4,:,2) = + dldHdz(:);

B(1,4,:,2) = - dldHdx(:);
B(2,3,:,2) = - dldHdx(:);

B(2,1,:,2) = - B(1,2,:,2);
B(4,3,:,2) = - B(3,4,:,2);

B(4,1,:,2) = - B(1,4,:,2);
B(3,2,:,2) = - B(2,3,:,2);


% m = 3
B(1,2,:,3) = - dldHdy(:);
B(3,4,:,3) = - dldHdy(:);

B(1,3,:,3) = + dldHdx(:);
B(2,4,:,3) = - dldHdx(:);

B(2,1,:,3) = - B(1,2,:,3);
B(4,3,:,3) = - B(3,4,:,3);

B(3,1,:,3) = - B(1,3,:,3);
B(4,2,:,3) = - B(2,4,:,3);



% factor of 1/6
B = 1./6. * B;

end

