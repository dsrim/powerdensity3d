function computeHij(sigma_type,const,N)
%COMPUTEHIJ
%
% load conductivity solutions then compute the power densities Hij
% requires outputs from computeIsoSols.
%
% saves output filename Hij_[sigma_type]_[const]_[N].mat

sol_fname = getSaveFname('IsoSols',sigma_type,const,N);
load(sol_fname);
varname = ['x', 'y', 'z'];
for i = 1:3
    for j = 1:3
        eval(['H' num2str(i) num2str(j) ' = zeros(N,N,N);'])
        eval(['H' num2str(i) num2str(j) ' = ' ...
                  'sigval.*(du' num2str(i) 'dx.*du' num2str(j) 'dx' ...
                         '+ du' num2str(i) 'dy.*du' num2str(j) 'dy' ...
                         '+ du' num2str(i) 'dz.*du' num2str(j) 'dz);'])
        eval(['S' num2str(i) num2str(j) ' = zeros(N,N,N);'])
        eval(['S' num2str(i) num2str(j) ' = ' ...
            'sqrt(sigval).*(du' num2str(i) 'd' varname(j) ');'])        
    end
end

% save to output file
save_fname = getSaveFname('Hij',sigma_type,const,N);
save(save_fname, ...
      'H11','H12','H13','H21','H22','H23','H31','H32','H33', ...
      'S11','S12','S13','S21','S22','S23','S31','S32','S33','h','-v7.3');


end
