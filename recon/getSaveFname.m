function save_fname = getSaveFname(name,sigma_type,const,N)
%GEATSAVEFNAME
%
% load the file named
%       [name]_[sigma_type]_[const]_[N].mat
% in the directory output_dir
%
% name, sigma_type: strings
% const, N:         numbers

output_dir = '../output/';
save_fname = [output_dir name '_' sigma_type ...
                              '_' num2str(const) '_' num2str(N) '.mat'];

end