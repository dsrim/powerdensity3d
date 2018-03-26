%% Run forward solver and reconstruction for anisotropic conductivities

addpath(genpath('../'));        % assume run directory is run/

N = 128;         % set grid size
hmax = 0.05;     % maximum mesh width
%const = 10;     % intensity of the anisotropy

gamma_list = {'anisophi','isotori','anisotori'};
          
% main loop 
for j = 1:length(gamma_list)
    for const = [2,20]
    gamma_type = gamma_list{j};     % select anisotropic conductivity gamma
        recon3p2(gamma_type,N,const,hmax)        % 3+2 reconstruction algo
        reconStab3p2(gamma_type,N,const,hmax)    % stabilized 3+2 recon. algo.
    end
end

