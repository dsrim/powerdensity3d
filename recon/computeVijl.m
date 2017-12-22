function computeVijl(sigma_type,const,N,varargin)
%computeVijl
%
% computes V_{ij}^l := \partial_l (t_{ik}) t^{kj} where T^T T = H^{-1}
%

Tij_fname  = getSaveFname('Tij',sigma_type,const,N);

disp(['(' mfilename ') loading Tij ...'])
load(Tij_fname);

opt = 'meshgrid';     % meshgrid option text
h = 2/(N-1);
for i = 1:3
  for j = 1:3
    for l = 1:3
      eval(['V' num2str(i) num2str(j) num2str(l) ' = zeros(N,N,N);'])
    end
    strij = [num2str(i) num2str(j)];
    eval(['[dT' strij 'dx1, dT' strij 'dx2, dT' strij 'dx3] = computeGrad(T' strij ',h,opt);'])
  end
end


for k = 1:(N^3)
  
  A  = [ dT11dx1(k),dT12dx1(k),dT13dx1(k);
         dT21dx1(k),dT22dx1(k),dT23dx1(k);
         dT31dx1(k),dT32dx1(k),dT33dx1(k) ];
  B  = [ t11(k),t12(k),t13(k);
         t21(k),t22(k),t23(k);
         t31(k),t32(k),t33(k) ];
  V = A*B;
  V111(k) = V(1,1);  V121(k) = V(1,2);  V131(k) = V(1,3);
  V211(k) = V(2,1);  V221(k) = V(2,2);  V231(k) = V(2,3);
  V311(k) = V(3,1);  V321(k) = V(3,2);  V331(k) = V(3,3);

  A  = [ dT11dx2(k),dT12dx2(k),dT13dx2(k);
         dT21dx2(k),dT22dx2(k),dT23dx2(k);
         dT31dx2(k),dT32dx2(k),dT33dx2(k) ];
  B  = [ t11(k),t12(k),t13(k);
         t21(k),t22(k),t23(k);
         t31(k),t32(k),t33(k) ];
  V = A*B;
  V112(k) = V(1,1);  V122(k) = V(1,2);   V132(k) = V(1,3);
  V212(k) = V(2,1);  V222(k) = V(2,2);   V232(k) = V(2,3);
  V312(k) = V(3,1);  V322(k) = V(3,2);   V332(k) = V(3,3);

  A  = [ dT11dx3(k),dT12dx3(k),dT13dx3(k);
         dT21dx3(k),dT22dx3(k),dT23dx3(k);
         dT31dx3(k),dT32dx3(k),dT33dx3(k) ];
  B  = [ t11(k),t12(k),t13(k);
         t21(k),t22(k),t23(k);
         t31(k),t32(k),t33(k) ];
  V = A*B;
  
  V113(k) = V(1,1);  V123(k) = V(1,2);  V133(k) = V(1,3);
  V213(k) = V(2,1);  V223(k) = V(2,2);  V233(k) = V(2,3);
  V313(k) = V(3,1);  V323(k) = V(3,2);  V333(k) = V(3,3);
end

%% Save variables

Vijl_fname = getSaveFname('Vijl',sigma_type,const,N);
disp(['(' mfilename ') saving to Vijl ...'])
save(Vijl_fname, ...
  'V111','V121','V131','V112','V122','V132','V113','V123','V133',...
  'V211','V221','V231','V212','V222','V232','V213','V223','V233',...
  'V311','V321','V331','V312','V322','V332','V313','V323','V333');




end

