function V = getVijl(sigma_name,const,N)
%GETVIJL
%
% compute / load Vijl then return it in array form.
% essentially a wrapper for computeVijl

Vijl_fname = getSaveFname('Vijl',sigma_name,const,N);
load(Vijl_fname)

M = N^3;
V = zeros(3,3,3,M);

% lazy fix to avoid nested MATLAB loops
% T is cholesky factor, lower diag *is admittedly* redundant
for k = 1:M
    
   V(1,1,1,k) = V111(k);
   V(1,1,2,k) = V112(k);
   V(1,1,3,k) = V113(k);
   V(1,2,1,k) = V121(k);
   V(1,2,2,k) = V122(k);
   V(1,2,3,k) = V123(k);
   V(1,3,1,k) = V131(k);
   V(1,3,2,k) = V132(k);
   V(1,3,3,k) = V133(k);
                 
   V(2,1,1,k) = V211(k);
   V(2,1,2,k) = V212(k);
   V(2,1,3,k) = V213(k);
   V(2,2,1,k) = V221(k);
   V(2,2,2,k) = V222(k);
   V(2,2,3,k) = V223(k);
   V(2,3,1,k) = V231(k);
   V(2,3,2,k) = V232(k);
   V(2,3,3,k) = V233(k);
                
   V(3,1,1,k) = V311(k);
   V(3,1,2,k) = V312(k);
   V(3,1,3,k) = V313(k);
   V(3,2,1,k) = V321(k);
   V(3,2,2,k) = V322(k);
   V(3,2,3,k) = V323(k);
   V(3,3,1,k) = V331(k);
   V(3,3,2,k) = V332(k);
   V(3,3,3,k) = V333(k);


end


end
