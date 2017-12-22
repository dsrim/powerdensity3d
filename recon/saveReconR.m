function saveReconR(sigma_type,const,N,R,error4R)
%SAVERECONR
%
% helper function to save R of size 3 x 3 x N^3
% to individual arrays reconRij of size N x N x N
%

reconR11 = reshape(R(1,1,:),N,N,N);
reconR12 = reshape(R(1,2,:),N,N,N);
reconR13 = reshape(R(1,3,:),N,N,N);

reconR21 = reshape(R(2,1,:),N,N,N);
reconR22 = reshape(R(2,2,:),N,N,N);
reconR23 = reshape(R(2,3,:),N,N,N);

reconR31 = reshape(R(3,1,:),N,N,N);
reconR32 = reshape(R(3,2,:),N,N,N);
reconR33 = reshape(R(3,3,:),N,N,N);

save_fname = getSaveFname('reconR',sigma_type,const,N);

save(save_fname, 'reconR11','reconR12','reconR13', ...
                 'reconR21','reconR22','reconR23', ...
                 'reconR31','reconR32','reconR33','error4R')
end
