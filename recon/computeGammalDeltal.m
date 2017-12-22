function computeGammalDeltal(sigma_name,resolution)
% 1 (1,2)
% 2 (2,3)
% 3 (1,3)


load(['LambdalFl_' sigma_name '_' num2str(resolution) '.mat'])
load(['Qijl_' sigma_name '_' num2str(resolution) '.mat'])
% load(['Tij_' sigma_name '_' num2str(resolution) '.mat']
load(['dQijldx_' sigma_name '_' num2str(resolution) '.mat'])

N = resolution;
Gm1 = zeros(N,N,N);
Gm2 = zeros(N,N,N);
Gm3 = zeros(N,N,N);
D11 = zeros(N,N,N);
D21 = zeros(N,N,N);
D31 = zeros(N,N,N);
D12 = zeros(N,N,N);
D22 = zeros(N,N,N);
D32 = zeros(N,N,N);
D13 = zeros(N,N,N);
D23 = zeros(N,N,N);
D33 = zeros(N,N,N);

for k = 1:N^3
  %% l = 1
  Q  = [Q111(k), Q121(k), Q131(k);
        Q211(k), Q221(k), Q231(k);
        Q311(k), Q321(k), Q331(k)];
 
  dQ = [dQdx111(k), dQdx121(k), dQdx131(k);
        dQdx211(k), dQdx221(k), dQdx231(k);
        dQdx311(k), dQdx321(k), dQdx331(k)];
      
  F = [ 0,          F1(k), F3(k);
       -F1(k),  0,         F2(k);
       -F3(k), -F2(k), 0        ];
     
  sksym = dQ*Q' + Q*F*Q';
  
  sksymPart = .5 * (sksym - sksym');   
  symPart = .5 * (sksym + sksym');   
  
  Gm1(k) = sksymPart(1,2);
  Gm2(k) = sksymPart(2,3);
  Gm3(k) = sksymPart(1,3);
  
  D11(k) = symPart(1,2);
  D21(k) = symPart(2,3);
  D31(k) = symPart(1,3);

  %% l = 2
  Q = [Q112(k), Q122(k), Q132(k);
       Q212(k), Q222(k), Q232(k);
       Q312(k), Q322(k), Q332(k)];
     
  dQ = [dQdx112(k), dQdx122(k), dQdx132(k);
        dQdx212(k), dQdx222(k), dQdx232(k);
        dQdx312(k), dQdx322(k), dQdx332(k)];
      
  sksym = dQ*Q'; 
  symPart = .5 * (sksym + sksym');   
  
  D12(k) = symPart(1,2);
  D22(k) = symPart(2,3);
  D32(k) = symPart(1,3);

  %% l = 3
  Q  = [Q113(k), Q123(k), Q133(k);
        Q213(k), Q223(k), Q233(k);
        Q313(k), Q323(k), Q333(k)];
      
  dQ = [dQdx113(k), dQdx123(k), dQdx133(k);
        dQdx213(k), dQdx223(k), dQdx233(k);
        dQdx313(k), dQdx323(k), dQdx333(k)];
      
  sksym = dQ*Q';
  symPart = .5 * (sksym + sksym');  
  
  D13(k) = symPart(1,2);
  D23(k) = symPart(2,3);
  D33(k) = symPart(1,3);

end
    
%% Save variables

output_dir = '../output/';
display(['(' mfilename ') plotting / saving outputs'])

save([output_dir 'GammalDeltal_' sigma_name '_' num2str(resolution) '.mat'], ...
              'Gm1','Gm2','Gm3',...
              'D11','D12','D13',... 
              'D21','D22','D23',... 
              'D31','D32','D33');

end
