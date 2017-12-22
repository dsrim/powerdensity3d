function plotEV(pos,components,sigmatype,N,m0,type)

if ((length(type) == 5) && prod(type == 'recon'))
    load(['reconRQ_' sigmatype '_' num2str(N) '.mat'])
    
    R11 = reconR11;
    R12 = reconR12;
    R13 = reconR13;
    
    R21 = reconR21;
    R22 = reconR22;
    R23 = reconR23;
    
    R31 = reconR31;
    R32 = reconR32;
    R33 = reconR33;
    
else
    load(['Tij_' sigmatype '_' num2str(N) '.mat'])
end
load(['sol_' sigmatype '_' num2str(N) '.mat'],'X','Y','Z')


s = 0.2;
EV = zeros(3,3,N);
loc = zeros(3,N);
i = pos(1); j = pos(2);
colors = ['r','b','k'];
for k = 0:(N-1)
% meshgrid ordering
  EV(:,:,k+1) = [R11(i+j*N^2+k*N) R12(i+j*N^2+k*N) R13(i+j*N^2+k*N);
                 R21(i+j*N^2+k*N) R22(i+j*N^2+k*N) R23(i+j*N^2+k*N);
                 R31(i+j*N^2+k*N) R32(i+j*N^2+k*N) R33(i+j*N^2+k*N)];
  loc(:,k+1) = [X(i+j*N^2+k*N) Y(i+j*N^2+k*N) Z(i+j*N^2+k*N)];
end

I0 = eye(3);

for m = components 
    
    ev1 = reshape(EV(1,m,:),1,N);
    ev2 = reshape(EV(2,m,:),1,N);
    ev3 = reshape(EV(3,m,:),1,N);    
    em = I0(:,m);
    for j = 1:m0:N
        q = quiver3(loc(1,j),loc(2,j),loc(3,j),ev1(j),ev2(j),ev3(j));
        aa = (1-min(4.*norm([ev1(j);ev2(j);ev3(j)] - em),1.));
        q.Color = aa*(.9*[1,1,1]) + (1-aa)*em';
        q.ShowArrowHead = 'off';
        q.AutoScaleFactor = 0.05;
        q.LineWidth = 2;
        hold on;  
        drawnow;
    end
end
% hold off;
axis equal
title(['eigenvectors, N= ' num2str(N) ', (i,j)=(' num2str(pos(1)) ', ' num2str(pos(2)) ')'])
xlabel('x'); ylabel('y'); zlabel('z');

end
