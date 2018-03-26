function computeSigmaPoisson(sigma_type,const,N)
%COMPUTESIGMAPOISSON
%
% solve a Poisson problem 

% load required data
disp(['(' mfilename ') loading Hij/reconR/Vijl ...']);
Hij_fname    = getSaveFname('Hij',   sigma_type,const,N);
Vijl_fname   = getSaveFname('Vijl',  sigma_type,const,N);
reconR_fname = getSaveFname('reconR',sigma_type,const,N);

load(Hij_fname)
load(reconR_fname)
load(Vijl_fname)


% compute source term (to the Poisson problem)
detH = zeros(N,N,N);

for k = 1:N^3
    H = [H11(k) H12(k) H13(k);
         H21(k) H22(k) H23(k);
         H31(k) H32(k) H33(k)];
    detH(k) = det(H);
end
[D1logdetH, D2logdetH, D3logdetH] = computeGrad(log(detH),h,'meshgrid');

VR1 = zeros(N,N,N);
VR2 = zeros(N,N,N);
VR3 = zeros(N,N,N);

for k = 1:N^3
  VRmat = zeros(3,1);
  symV = zeros(3,3,3);

  rR = [reconR11(k), reconR12(k), reconR13(k);
        reconR21(k), reconR22(k), reconR23(k);
        reconR31(k), reconR32(k), reconR33(k)];

  V1 = [V111(k) V121(k) V131(k); 
        V211(k) V221(k) V231(k); 
        V311(k) V321(k) V331(k)]; % x-component of V
  V2 = [V112(k) V122(k) V132(k); 
        V212(k) V222(k) V232(k); 
        V312(k) V322(k) V332(k)]; % y-component of V
  V3 = [V113(k) V123(k) V133(k); 
        V213(k) V223(k) V233(k); 
        V313(k) V323(k) V333(k)]; % z-component of V

  symV(:,1,:) =  .5*(V1 + V1');
  symV(:,2,:) =  .5*(V2 + V2');
  symV(:,3,:) =  .5*(V3 + V3');
  
  VRmat = [trace([symV(1,:,1); symV(1,:,2); symV(1,:,3)]*rR);
           trace([symV(2,:,1); symV(2,:,2); symV(2,:,3)]*rR);
           trace([symV(3,:,1); symV(3,:,2); symV(3,:,3)]*rR)];

  VR1(k) = rR(1,:)*VRmat;
  VR2(k) = rR(2,:)*VRmat;
  VR3(k) = rR(3,:)*VRmat;
end

[dVRdx,     ~,     ~] = computeGrad(2.*VR1 + .5*D1logdetH,h,'meshgrid');
[    ~, dVRdy,     ~] = computeGrad(2.*VR2 + .5*D2logdetH,h,'meshgrid');
[    ~,     ~, dVRdz] = computeGrad(2.*VR3 + .5*D3logdetH,h,'meshgrid');

b = 2/3*(dVRdx + dVRdy + dVRdz);    % compute RHS


% solve Poisson Problem
u = zeros(N,N,N);
J = getBdryAll(N); J = J(:,1);
logsigmatype = ['log' sigma_type];
w0 = -1; w1 =1;
x = linspace(w0,w1,N);
[X,Y,Z] = meshgrid(x,x,x);

u(J)     = iso_conductivity(logsigmatype,const,[X(J),Y(J),Z(J)]');
logsigma = solveScalarConductivity(u,b,'identity',const);
sigma    = exp(logsigma);

% compute error
sol_fname = getSaveFname('IsoSols',sigma_type,const,N);
load(sol_fname, 'sigval');
e4sigma = zeros(N,N,N);
e4sigma = sigval(:) - sigma(:);

% Save variables
disp(['(' mfilename ') saving outputs ']);
save_fname = getSaveFname('reconSigma',sigma_type,const,N);
save(save_fname,'logsigma','sigma','e4sigma','-v7.3');

end






