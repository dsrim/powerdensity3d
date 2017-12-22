function error4Sigma = computeError4Sigma(sigma_type,const,N)

% load data
recon_fname = getSaveFname('reconSigma',sigma_type,const,N);
sol_fname   = getSaveFname('IsoSols',   sigma_type,const,N);

disp(['(' mfilename ') loading reconR ...'])
load(recon_fname,'sigma');
load(sol_fname,'sigval');

error4Sigma = zeros(N,N,N);
error4Sigma(:) = sigma(:) - sigval(:);

end
