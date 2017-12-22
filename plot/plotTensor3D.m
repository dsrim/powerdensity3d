function plotTensor3D(T,m,varargin)
%plotTensor3D plots 3D slice of a 3x3 tensor.
%
% input T should be a 3 x 3 x M array where M is the # of 3D gridpts
% plotTensor3D(T,m,varargin)
% to get a "upser" title, requires bioinformatics toolbox routine suptitle
%
%

if ~isempty(varargin)
    varname = varargin{1};
    if length(varargin) > 1
    plotname = varargin{2};
    end
else
    varname = 'T';
end

M = max(size(T));
N = floor(M^(1/3));
if (M ~= N^3)
  N = N + 1;
end

vmin = min(T(:));
vmax = max(T(:));
subplot(331); ezplot3D(reshape(T(1,1,:),N,N,N),m); title([varname '_1_1'])
caxis([vmin,vmax]);
subplot(332); ezplot3D(reshape(T(1,2,:),N,N,N),m); title([varname '_1_2'])
caxis([vmin,vmax]);
subplot(333); ezplot3D(reshape(T(1,3,:),N,N,N),m); title([varname '_1_3'])
caxis([vmin,vmax]);

subplot(334); ezplot3D(reshape(T(2,1,:),N,N,N),m); title([varname '_2_1'])
caxis([vmin,vmax]);
subplot(335); ezplot3D(reshape(T(2,2,:),N,N,N),m); title([varname '_2_2'])
caxis([vmin,vmax]);
subplot(336); ezplot3D(reshape(T(2,3,:),N,N,N),m); title([varname '_2_3'])
caxis([vmin,vmax]);

subplot(337); ezplot3D(reshape(T(3,1,:),N,N,N),m); title([varname '_3_1'])
caxis([vmin,vmax]);
subplot(338); ezplot3D(reshape(T(3,2,:),N,N,N),m); title([varname '_3_2'])
caxis([vmin,vmax]);
subplot(339); ezplot3D(reshape(T(3,3,:),N,N,N),m); title([varname '_3_3'])
caxis([vmin,vmax]);

if exist('plotname','var') && (exist('suptitle','file') ~= 0)
    suptitle(plotname);
end

end
