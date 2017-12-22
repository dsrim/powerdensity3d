function reconIso(sigma_type, const, N)
%RECONISO
%
% run reconstruction for scalar conductivity
% by solving a dynamical system for unknown SO(3) using quaternions,
% then by solving a poisson problem that follows.

h = 2/(N-1);                    % grid width

%% generate three solutions to the conductivity problem

sol_fname = getSaveFname('IsoSols',sigma_type,const,N);
if exist(sol_fname,'file')
    disp(['saved file found: ' sol_fname])
else
    computeIsoSols(sigma_type,const,N);
end

%% compute power density Hij 

Hij_fname = getSaveFname('Hij',sigma_type,const,N);
if exist(Hij_fname, 'file')
    disp(['saved file found: ' Hij_fname])
else
    computeHij(sigma_type,const,N);
end

%% compute cholesky factor Tij 

Tij_fname = getSaveFname('Tij',sigma_type,const,N);
if exist(Tij_fname, 'file')
    disp(['saved file found: ' Tij_fname])
else
    computeTij(sigma_type,const,N);
end

%% compute Vijl

Vijl_fname = getSaveFname('Vijl',sigma_type,const,N);
if exist(Vijl_fname, 'file')
    disp(['saved file found: ' Vijl_fname])
else
    computeVijl(sigma_type,const,N);
end


%% solve dynamical system in quaternions

X_fname = getSaveFname('X',sigma_type,const,N);
if exist(X_fname, 'file')
    disp(['saved file found: ' X_fname])
else
    computeXall(sigma_type,const,N);
end

%% check errors
load(X_fname)
R = computeR3(Xall);        % convert to 3x3 rotation matrix
Rtrue = getRtrue(sigma_type,const,N);
error4R = norm(R(:) - Rtrue(:))*h^(3/2);
display(['L2 error for R is: ' num2str(error4R) ])

%% save R if not already saved
R_fname = getSaveFname('reconR',sigma_type,const,N);
if exist(R_fname, 'file')
    disp(['saved file found: ' R_fname])
else
    saveReconR(sigma_type,const,N,R,error4R);
end

%% solve Poisson problem to solve for log(sigma)
reconSigma_fname = getSaveFname('reconSigma',sigma_type,const,N);
if exist(reconSigma_fname, 'file')
    disp(['saved file found: ' reconSigma_fname])
else
    computeSigmaPoisson(sigma_type,const,N);
end

%% compute error
error4sigma = computeError4Sigma(sigma_type,const,N);


