function plotdEV(pos,components,sigmatype,N)

load(['dQijdx_' sigmatype '_' num2str(N) '.mat'])
load(['sol_' sigmatype '_' num2str(N) '.mat'],'X','Y','Z')

EV = zeros(3,3,N);
loc = zeros(3,N);
i = pos(1); j = pos(2);
colors = ['r','g','b'];
for k = 0:(N-1)
  EV(:,:,k+1) = [dQdx11(i+j*N^2+k*N) dQdx12(i+j*N^2+k*N) dQdx13(i+j*N^2+k*N);
                 dQdx21(i+j*N^2+k*N) dQdx22(i+j*N^2+k*N) dQdx23(i+j*N^2+k*N);
                 dQdx31(i+j*N^2+k*N) dQdx32(i+j*N^2+k*N) dQdx33(i+j*N^2+k*N)];
  loc(:,k+1) = [X(i+j*N^2+k*N) Y(i+j*N^2+k*N) Z(i+j*N^2+k*N)];
end

figure;
for m = components 
    ev1 = reshape(EV(1,m,:),1,N);
    ev2 = reshape(EV(2,m,:),1,N);
    ev3 = reshape(EV(3,m,:),1,N);
    quiver3(loc(1,:),loc(2,:),loc(3,:),ev1,ev2,ev3,colors(m));
    hold on;                                                       
end
hold off;
axis equal
title(['eigenvectors, N= ' num2str(N) ', (i,j)=(' num2str(pos(1)) ', ' num2str(pos(2)) ')'])
xlabel('x'); ylabel('y'); zlabel('z');

end
