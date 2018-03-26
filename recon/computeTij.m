function computeTij(sigma_type,const,N)
%COMPUTETIJ
%
% compute Cholesky factor Tij of Hij inverse
%

Hij_fname = getSaveFname('Hij',sigma_type,const,N);
load(Hij_fname);

for i = 1:3
    for j = 1:3
        eval(['T' num2str(i) num2str(j) ' = zeros(N,N,N);'])
        eval(['t' num2str(i) num2str(j) ' = zeros(N,N,N);'])
        eval(['R' num2str(i) num2str(j) ' = zeros(N,N,N);'])
    end
end

for k = 1:(N^3)

    H = [H11(k) H12(k) H13(k); 
         H21(k) H22(k) H23(k); 
         H31(k) H32(k) H33(k)];

    
    T = chol(H\eye(3));
    T11(k) = T(1,1);
    T12(k) = T(1,2);
    T13(k) = T(1,3);
    T21(k) = T(2,1);
    T22(k) = T(2,2);
    T23(k) = T(2,3);
    T31(k) = T(3,1);
    T32(k) = T(3,2);
    T33(k) = T(3,3);
    
    t = T \ eye(3);
    t11(k) = t(1,1);
    t12(k) = t(1,2);
    t13(k) = t(1,3);
    t21(k) = t(2,1);
    t22(k) = t(2,2);
    t23(k) = t(2,3);
    t31(k) = t(3,1);
    t32(k) = t(3,2);
    t33(k) = t(3,3);
    
    S = [S11(k) S12(k) S13(k); 
         S21(k) S22(k) S23(k); 
         S31(k) S32(k) S33(k)];   

    R = S*T';

    R11(k) = R(1,1);
    R12(k) = R(1,2);
    R13(k) = R(1,3);
    R21(k) = R(2,1);
    R22(k) = R(2,2);
    R23(k) = R(2,3);
    R31(k) = R(3,1);
    R32(k) = R(3,2);
    R33(k) = R(3,3);
end

%% Save variables
Tij_fname = getSaveFname('Tij',sigma_type,const,N);
save(Tij_fname, ...
      'T11','T12','T13','T21','T22','T23','T31','T32','T33',...
      't11','t12','t13','t21','t22','t23','t31','t32','t33',...
      'R11','R12','R13','R21','R22','R23','R31','R32','R33','h','-v7.3');

end
