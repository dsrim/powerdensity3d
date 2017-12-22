%% run options, grid-size and sigma

addpath(genpath('../'))

N = 64;
const = 1;
sigma_list = {'gaussians', 'smoothtori', 'balls','tinyballs','tori'};

for j = 1:length(sigma_list)

    sigma_type = sigma_list{j};
    reconIso(sigma_type, const, N);

end
