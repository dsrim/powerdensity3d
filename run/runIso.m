%% run options, grid-size and sigma

addpath(genpath('../'))

N = 64;
hmax = 0.1;
const = 2;
%sigma_list = {'one_gaussian','gaussians', 'smoothtori', 'balls','tinyballs','tori','smoothtori_exp1','smoothtori_exp1_alt'};
sigma_list = {'smoothtori_exp1'};

for j = 1:length(sigma_list)

    sigma_type = sigma_list{j};
    reconIso(sigma_type, const, N, hmax);

end
