function plotEigval(pos,components,sigmatype,N)

load(['LambdalFl_' sigmatype '_' num2str(N) '.mat'])
load(['sol_' sigmatype '_' num2str(N) '.mat'],'X','Y','Z')

EV = zeros(3,3,N);
loc = zeros(3,N);
i = pos(1); j = pos(2);
for k = 0:(N-1)
  EV(:,:,k+1) = [Lambda11(i+j*N^2+k*N) Lambda12(i+j*N^2+k*N) Lambda13(i+j*N^2+k*N);
                 Lambda21(i+j*N^2+k*N) Lambda22(i+j*N^2+k*N) Lambda23(i+j*N^2+k*N);
                 Lambda31(i+j*N^2+k*N) Lambda32(i+j*N^2+k*N) Lambda33(i+j*N^2+k*N)];
  loc(:,k+1) = [X(i+j*N^2+k*N) Y(i+j*N^2+k*N) Z(i+j*N^2+k*N)];
end

figure;clf;
% plot3(loc(1,:),loc(2,:),loc(3,:),'k');
% hold on;  
for m = components 
    ev1 = reshape(EV(1,m,:),1,N)
    ev2 = reshape(EV(2,m,:),1,N)
    ev3 = reshape(EV(3,m,:),1,N)
    plot(loc(1,:),ev1,'r-x',loc(1,:),ev2,'g-+',loc(1,:),ev3,'b-o'); hold on;
%     plot3(loc(1,:)+ev1,loc(2,:),loc(3,:),'r');
%     plot3(loc(1,:),loc(2,:)+ev2,loc(3,:),'g');
%     plot3(loc(1,:),loc(2,:),loc(3,:)+ev3,'b');                                                   
end
hold off;
axis equal
title(['eigenvalues, N= ' num2str(N) ', (i,j)=(' num2str(pos(1)) ', ' num2str(pos(2)) ')'])
xlabel('x'); ylabel('y'); zlabel('z');

end
