function computeSigma(sigmatype,N,l, varargin)



if ~isempty(varargin) && strcmp(varargin{1}, 'QT')
    display(['(' mfilename ') loading Hij/QreconR/QVijl ...'])
    load(['Hij_' sigmatype '_' num2str(N) '.mat'])
    load(['QreconR_' sigmatype '_' num2str(N) '.mat'])
    load(['QVijl_' sigmatype '_' num2str(N) '.mat'])   
else
    display(['(' mfilename ') loading Hij/reconR/Vijl ...'])
    load(['Hij_' sigmatype '_' num2str(N) '.mat'])
    load(['reconR_' sigmatype '_' num2str(N) '.mat'])
    load(['Vijl_' sigmatype '_' num2str(N) '.mat'])
end


detH = zeros(N,N,N);
for k = 1:N^3
    H = [H11(k) H12(k) H13(k);
         H21(k) H22(k) H23(k);
         H31(k) H32(k) H33(k)];
    detH(k) = det(H);
end

[D1logdetH, D2logdetH, D3logdetH] = computeGrad(log(detH),h);

h = 2/(N-1);
logsigma = zeros(N,N,N);

load(['sol_' sigmatype '_' num2str(N) '.mat'], 'sigval')

if l == 1
for i = 1:N
    for k = 1:N
        logsigma(i,1,k) = log(sigval(i,1,k));
    end
end


%load(['Tij_' sigmatype '_' num2str(N) '.mat'],...
%      'R11','R12','R13',...
%      'R21','R22','R23',...
%      'R31','R32','R33')


for i = 1:N
for j = 1:(N-1)
for k = 1:N
    R1 = [reconR11(i,j,k); reconR21(i,j,k); reconR31(i,j,k)];
    R2 = [reconR12(i,j,k); reconR22(i,j,k); reconR32(i,j,k)];
    R3 = [reconR13(i,j,k); reconR23(i,j,k); reconR33(i,j,k)];
logsigma(i,j+1,k) = logsigma(i,j,k) ...
+ h*2/3*(1/2*D1logdetH(i,j,k)   ...
    + 2*[V111(i,j,k), V112(i,j,k), V113(i,j,k)]*R1*R1(1)  ...
    + 2*[V221(i,j,k), V222(i,j,k), V223(i,j,k)]*R2*R2(1)  ...
    + 2*[V331(i,j,k), V332(i,j,k), V333(i,j,k)]*R3*R3(1)  ...
    + ([V211(i,j,k), V212(i,j,k), V213(i,j,k)] + [V121(i,j,k), V122(i,j,k), V123(i,j,k)])*R1*R2(1) ...
    + ([V211(i,j,k), V212(i,j,k), V213(i,j,k)] + [V121(i,j,k), V122(i,j,k), V123(i,j,k)])*R2*R1(1) ...
    + ([V311(i,j,k), V312(i,j,k), V313(i,j,k)] + [V131(i,j,k), V132(i,j,k), V133(i,j,k)])*R3*R1(1) ...
    + ([V311(i,j,k), V312(i,j,k), V313(i,j,k)] + [V131(i,j,k), V132(i,j,k), V133(i,j,k)])*R1*R3(1) ...
    + ([V321(i,j,k), V322(i,j,k), V323(i,j,k)] + [V231(i,j,k), V232(i,j,k), V233(i,j,k)])*R3*R2(1) ...
    + ([V321(i,j,k), V322(i,j,k), V323(i,j,k)] + [V231(i,j,k), V232(i,j,k), V233(i,j,k)])*R2*R3(1));
%   logsigma(i,j+1,k) = logsigma(i,j,k) + h*2/3*(1/2*D1detH(i,j,k));
end
end
end
elseif l == 2
  
  for i = 1:N
    for k = 1:N
        logsigma(1,i,k) = log(sigval(1,i,k));
    end
end


%load(['Tij_' sigmatype '_' num2str(N) '.mat'],...
%      'R11','R12','R13',...
%      'R21','R22','R23',...
%      'R31','R32','R33')


for i = 1:(N-1)
for j = 1:N
for k = 1:N
    R1 = [reconR11(i,j,k); reconR21(i,j,k); reconR31(i,j,k)];
    R2 = [reconR12(i,j,k); reconR22(i,j,k); reconR32(i,j,k)];
    R3 = [reconR13(i,j,k); reconR23(i,j,k); reconR33(i,j,k)];
logsigma(i+1,j,k) = logsigma(i,j,k) ...
+ h*2/3*(1/2*D1logdetH(i,j,k)   ...
    + 2*[V111(i,j,k), V112(i,j,k), V113(i,j,k)]*R1*R1(2)  ...
    + 2*[V221(i,j,k), V222(i,j,k), V223(i,j,k)]*R2*R2(2)  ...
    + 2*[V331(i,j,k), V332(i,j,k), V333(i,j,k)]*R3*R3(2)  ...
    + ([V211(i,j,k), V212(i,j,k), V213(i,j,k)] + [V121(i,j,k), V122(i,j,k), V123(i,j,k)])*R1*R1(2) ...
    + ([V311(i,j,k), V312(i,j,k), V313(i,j,k)] + [V131(i,j,k), V132(i,j,k), V133(i,j,k)])*R3*R1(2) ...
    + ([V321(i,j,k), V322(i,j,k), V323(i,j,k)] + [V231(i,j,k), V232(i,j,k), V233(i,j,k)])*R3*R2(2));

end
end
end

end

sigma = exp(logsigma);

error4sigma = zeros(N,N,N);
for k = 1:N^3
    error4sigma(k) = abs(sigval(k) - sigma(k));
end

%% Save variables

output_dir = '../output/';
display(['(' mfilename ') plotting / saving outputs'])



if ~isempty(varargin)
if strcmp(varargin{1}, 'QT')
    display(['(' mfilename ') saving to QreconSigma ...'])
    save([output_dir 'QreconSigma_' sigmatype '_' num2str(N) '.mat'],'logsigma','sigma','error4sigma');
end
else
    display(['(' mfilename ') saving to reconSigma ...'])
    save([output_dir 'reconSigma_' sigmatype '_' num2str(N) '.mat'],'logsigma','sigma','error4sigma');
end

end
