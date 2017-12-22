function T = getTij(sigma_name,N,output_dir)
%GETRTRUE
%
%  helper function that loads the "true" R matrix

Tij_fname = [output_dir 'Tij_' sigma_name '_' num2str(N) '.mat'];

if exist(Tij_fname, 'file')
    display(['saved file found, loading: ' Tij_fname])
    load(Tij_fname)
else
    computeTij(sigma_name,N);
end

T = zeros(3,3,N^3);

T(1,1,:) = T11(:);
T(1,2,:) = T12(:);
T(1,3,:) = T13(:);
 
T(2,1,:) = T21(:);
T(2,2,:) = T22(:);
T(2,3,:) = T23(:);
 
T(3,1,:) = T31(:);
T(3,2,:) = T32(:);
T(3,3,:) = T33(:);

end

