function plotFDM3D(X,Y,Z,x,m,figtitle,filename)
a = min(min(min((X))));
b = max(max(max((X))));
N = size(X,1);
xslice = linspace(a,b,m); xslice = xslice(2:(end-1));
yslice = linspace(a,b,m); yslice = yslice(2:(end-1));
zslice = linspace(a,b,m); zslice = zslice(2:(end-1));
figure;
x = reshape(x,N,N,N);
f = slice(X,Y,Z,x,xslice,yslice,zslice);
alphamap('rampdown'); alphamap('increase',.1)
axis equal; colorbar; title(figtitle);
alpha('color'); alphamap('rampdown'); alphamap('increase',.1)
xlabel('x')
ylabel('y')
zlabel('z')
saveas(gcf, ['results/' filename '_' num2str(N) '.fig'], 'fig')

end