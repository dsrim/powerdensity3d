function [Z,ZHOm] = computeNormalMatStab(Mu,H,N)
%COMPUTENORMALMATSTAB
%
% compute Z and ZHOm from Mu and H (stabilized)
%
% if H is of size (K-3) x (K-3) x M where M = N^3 is the no. of 3D gridpts
% minimal number of K would be K=2
%
% returns 
%   Z       of size     3 x 3 x K x M
%   ZHOm    of size 3 x 3 x 3 x K x M

h = 2/(N-1);        % bad
K = size(Mu,2);     % let's assume that K=2 for the moment
M = size(Mu,3);
      
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];

Om1 = e2*e3' - e3*e2';
Om2 = e3*e1' - e1*e3';
Om3 = e1*e2' - e2*e1';

Z = zeros(3,3,K,M);

detH = det3(H)';
[ddetHdx,ddetHdy,ddetHdz] = computeGrad(detH,h);

for k = 1:K
    
    [dmu1dx,dmu1dy,dmu1dz] = computeGrad(squeeze(Mu(1,k,:)),h);
    [dmu2dx,dmu2dy,dmu2dz] = computeGrad(squeeze(Mu(2,k,:)),h);
    [dmu3dx,dmu3dy,dmu3dz] = computeGrad(squeeze(Mu(3,k,:)),h);
    

    Z(1,1,k,:) = dmu1dx.*detH - squeeze(Mu(1,k,:)).*ddetHdx;
    Z(2,1,k,:) = dmu1dy.*detH - squeeze(Mu(1,k,:)).*ddetHdy;
    Z(3,1,k,:) = dmu1dz.*detH - squeeze(Mu(1,k,:)).*ddetHdz;
    
    Z(1,2,k,:) = dmu2dx.*detH - squeeze(Mu(2,k,:)).*ddetHdx;
    Z(2,2,k,:) = dmu2dy.*detH - squeeze(Mu(2,k,:)).*ddetHdy;
    Z(3,2,k,:) = dmu2dz.*detH - squeeze(Mu(2,k,:)).*ddetHdz;
    
    Z(1,3,k,:) = dmu3dx.*detH - squeeze(Mu(3,k,:)).*ddetHdx;
    Z(2,3,k,:) = dmu3dy.*detH - squeeze(Mu(3,k,:)).*ddetHdy;
    Z(3,3,k,:) = dmu3dz.*detH - squeeze(Mu(3,k,:)).*ddetHdz;

end


ZHOm = zeros(3,3,3,K,M);

for k = 1:K
  
  ZHOm(:,:,1,k,:) = squeeze(...
                    matmul3(matmul3(Z(:,:,k,:), H(1:3,1:3,:)), Om1));

  ZHOm(:,:,2,k,:) = squeeze(...
                    matmul3(matmul3(Z(:,:,k,:), H(1:3,1:3,:)), Om2));

  ZHOm(:,:,3,k,:) = squeeze(...
                    matmul3(matmul3(Z(:,:,k,:), H(1:3,1:3,:)), Om3));
                 
                  
end


end
