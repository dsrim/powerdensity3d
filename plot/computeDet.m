for N = [40]
close all
%N = 100;
rsNo = N^2;
tol = 1e-8;
sigmatype = 'smoothtori';
w0 = -1; w1 = 1;
addpath('../reconst')

[h, x3d1, x3d2, x3d3, sigval, iX, iY, iZ, X, Y, Z] = FDM2Shepp3D(N,rsNo,tol,sigmatype);


% compute gradients (1st order)
du1dx = (x3d1(:,[2:end,end],:) - x3d1(:,[1:(end-1), (end-1)],:))/h;
du2dx = (x3d2(:,[2:end,end],:) - x3d2(:,[1:(end-1), (end-1)],:))/h;
du3dx = (x3d3(:,[2:end,end],:) - x3d3(:,[1:(end-1), (end-1)],:))/h;
du1dy = (x3d1([2:end,end],:,:) - x3d1([1:(end-1), (end-1)],:,:))/h;
du2dy = (x3d2([2:end,end],:,:) - x3d2([1:(end-1), (end-1)],:,:))/h;
du3dy = (x3d3([2:end,end],:,:) - x3d3([1:(end-1), (end-1)],:,:))/h;
du1dz = (x3d1(:,:,[2:end,end]) - x3d1(:,:,[1:(end-1), (end-1)]))/h;
du2dz = (x3d2(:,:,[2:end,end]) - x3d2(:,:,[1:(end-1), (end-1)]))/h;
du3dz = (x3d3(:,:,[2:end,end]) - x3d3(:,:,[1:(end-1), (end-1)]))/h;

detS = zeros(N^3,1);
detDU = zeros(N^3,1);

for k = 1:N^3
   S = sqrt(sigval(k)).*[du1dx(k), du1dy(k), du1dz(k);
                         du2dx(k), du2dy(k), du2dz(k);
                         du3dx(k), du3dy(k), du3dz(k)];
  DU = [du1dx(k), du1dy(k), du1dz(k);
      du2dx(k), du2dy(k), du2dz(k);
      du3dx(k), du3dy(k), du3dz(k)];                  
  detS(k,1) = det(S);
  detDU(k,1) = det(DU);
end

[mindet,minI] = min(detS);
[maxdet,maxI] = max(detS);
display(['min/max det(S): ' num2str(mindet) ' / ' num2str(maxdet)])


xslice = (w0 + .5):0.5:(w1 - .5);
yslice = (w0 + .5):0.5:(w1 - .5);
zslice = (w0 + .5):0.5:(w1 - .5);

figure;
detS3d = reshape(detS,N,N,N);
g = slice(X,Y,Z,detS3d,xslice,yslice,zslice);
alphamap('rampdown'); alphamap('increase',.1)
axis equal; colorbar; title('det(S)');
xlabel('x'); ylabel('y'); zlabel('z');
saveas(gcf, ['detS_' num2str(N) '.fig'], 'fig')

figure;
detDU3d = reshape(detDU,N,N,N);
g = slice(X,Y,Z,detDU3d,xslice,yslice,zslice);
alphamap('rampdown'); alphamap('increase',.1)
axis equal; colorbar; title('det(DU)');
xlabel('x'); ylabel('y'); zlabel('z');
saveas(gcf, ['detDU_' num2str(N) '.fig'], 'fig')



save(['data_' sigmatype '_' num2str(N) '.mat'],'x3d1','x3d2','x3d3',...
                                'du1dx','du2dx','du3dx',...
                                'du1dy','du2dy','du3dy',...
                                'du1dz','du2dz','du3dz',...
                                'detS','X','Y','Z','h','detDU','sigval')

end
