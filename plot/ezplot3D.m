function ezplot3D(W,M,varargin)
N = size(W,1);        % assuming equal mesh-size in x,y,z directions
w0 = -1; w1 = 1;

om = 3;
x = linspace(w0,w1,N); 
x = x((1+om):(end-om)); % omit endpts

[X,Y,Z] = meshgrid(x,x,x);  % slice() accepts only meshgrids
if ~isempty(varargin) && strcmp(varargin{1},'meshgrid')
else
    W = permute(W,[2,1,3]);
end

d=.5*3/M;

if M > 1
    xslice = linspace((w0 + d),(w1 - d),M);
    yslice = linspace((w0 + d),(w1 - d),M);
    zslice = linspace((w0 + d),(w1 - d),M);
else
    xslice = (w0 + w1)/2;
    yslice = (w0 + w1)/2;
    zslice = (w0 + w1)/2; 
end

g = slice(X,Y,Z,W((1+om):(end-om),(1+om):(end-om),(1+om):(end-om)),xslice,yslice,zslice); %omit endpts
xlabel('x'); ylabel('y'); zlabel('z');

% set(g,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
if N >= 80
  set(g,'EdgeColor','none')
end

% alpha('color'); alphamap('rampdown'); alphamap('increase',.2)
axis equal; colorbar; 
if nargin ==3
  title(varargin{1})
end

shading flat

end
