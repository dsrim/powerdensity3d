function [Z,ZHOm] = computeNormalMatPermute(Mu,H,N,sInd)
%COMPUTENORMALMATPERMUTE
%
% compute Z and ZHOm from Mu and H
%
% if H is of size (K-3) x (K-3) x M where M = N^3 is the no. of 3D gridpts
% minimal number of K would be K=2
%
% returns 
%   Z       of size     3 x 3 x K x M
%   ZHOm    of size 3 x 3 x 3 x K x M

h = 2/(N-1);        % bad
K = size(Mu,2);     % let's assume that K=2 for the moment
M = max(size(Mu));
      
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];

Om1 = e2*e3' - e3*e2';
Om2 = e3*e1' - e1*e3';
Om3 = e1*e2' - e2*e1';

Z = zeros(3,3,K,M);


for k = 1:K
    kstr = num2str(k);
    
    eval(['[dmu1dx,dmu1dy,dmu1dz] = computeGrad(Mu(1,' kstr ',:),h);'])
    eval(['[dmu2dx,dmu2dy,dmu2dz] = computeGrad(Mu(2,' kstr ',:),h);'])
    eval(['[dmu3dx,dmu3dy,dmu3dz] = computeGrad(Mu(3,' kstr ',:),h);'])
    
    eval(['Z(1,1,' kstr ',:) = dmu1dx;'])
    eval(['Z(2,1,' kstr ',:) = dmu1dy;'])
    eval(['Z(3,1,' kstr ',:) = dmu1dz;'])
    
    eval(['Z(1,2,' kstr ',:) = dmu2dx;'])
    eval(['Z(2,2,' kstr ',:) = dmu2dy;'])
    eval(['Z(3,2,' kstr ',:) = dmu2dz;'])
    
    eval(['Z(1,3,' kstr ',:) = dmu3dx;'])
    eval(['Z(2,3,' kstr ',:) = dmu3dy;'])
    eval(['Z(3,3,' kstr ',:) = dmu3dz;'])

end


ZHOm = zeros(3,3,3,K,M);

parfor j = 1:M
  sIndx = sInd(:,j);
  for k = 1:K
    Hx = H(sIndx(1:3),sIndx(1:3),j);
    Zx = Z(:,:,k,j)
    
    ZHOm(:,:,1,k,j) = Zx*Hx*Om1;
    ZHOm(:,:,2,k,j) = Zx*Hx*Om2;
    ZHOm(:,:,3,k,j) = Zx*Hx*Om3;

  end
end


end
