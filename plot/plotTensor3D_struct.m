function plotTensor3D_struct(T,m,varargin)
%plotTensor3D plots 3D slice of a 3x3 tensor.
%
% requires bioinformatics toolbox routine suptitle

if ~isempty(varargin)
    varname = varargin{1};
    if length(varargin) > 1
    plotname = varargin{2};
    end
else
    varname = 'T';
end

subplot(331); ezplot3D(T.m11,m); title([varname '_1_1'])
subplot(332); ezplot3D(T.m12,m); title([varname '_1_2'])
subplot(333); ezplot3D(T.m13,m); title([varname '_1_3'])

subplot(334); ezplot3D(T.m21,m); title([varname '_2_1'])
subplot(335); ezplot3D(T.m22,m); title([varname '_2_2'])
subplot(336); ezplot3D(T.m23,m); title([varname '_2_3'])

subplot(337); ezplot3D(T.m31,m); title([varname '_3_1'])
subplot(338); ezplot3D(T.m32,m); title([varname '_3_2'])
subplot(339); ezplot3D(T.m33,m); title([varname '_3_3'])

if exist('plotname','var')
    suptitle(plotname);
end

end
