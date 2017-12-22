function S = getSij(sigma_name,N,output_dir)
%GETRTRUE
%
%  helper function that loads the "true" R matrix

Hij_fname = [output_dir 'Hij_' sigma_name '_' num2str(N) '.mat']

if exist(Hij_fname, 'file')
    display(['saved file found, loading: ' Hij_fname])
    load(Hij_fname)
else
    computeHij(sigma_name,N);
end

S = zeros(3,3,N^3);

S(1,1,:) = S11(:);
S(1,2,:) = S12(:);
S(1,3,:) = S13(:);
 
S(2,1,:) = S21(:);
S(2,2,:) = S22(:);
S(2,3,:) = S23(:);
 
S(3,1,:) = S31(:);
S(3,2,:) = S32(:);
S(3,3,:) = S33(:);

end


