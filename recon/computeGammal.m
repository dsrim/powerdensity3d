function computeGammal(sigma_name,resolution)
% 1 (1,2)
% 2 (2,3)
% 3 (1,3)


load(['LambdalFl_' sigma_name '_' num2str(resolution) '.mat'])
load(['Qij_' sigma_name '_' num2str(resolution) '.mat'])
load(['dQijdx_' sigma_name '_' num2str(resolution) '.mat'])

N = resolution;
Gammam1 = zeros(N,N,N);
Gammam2 = zeros(N,N,N);
Gammam3 = zeros(N,N,N);
skError = zeros(N,N,N);
skErrorSym = zeros(N,N,N);
skErrorDiag = zeros(N,N,N);
I = eye(3);
h = 2/(N-1);
for k = 1:N^3
  Q = [Q11(k), Q12(k), Q13(k);
       Q21(k), Q22(k), Q23(k);
       Q31(k), Q32(k), Q33(k)];
  Qp1= Q;
  if k < (N^3 - N)
  Qp1 = [Q11(k+N), Q12(k+N), Q13(k+N);
         Q21(k+N), Q22(k+N), Q23(k+N);
         Q31(k+N), Q32(k+N), Q33(k+N)];
  end
  chg = norm(.5/h*(Qp1*Q' + Q*Qp1' - 2*I),'fro');
  dQ = [dQdx11(k), dQdx12(k), dQdx13(k);
        dQdx21(k), dQdx22(k), dQdx23(k);
        dQdx31(k), dQdx32(k), dQdx33(k)];
  F = [ 0,          F1(k), F3(k);
       -F1(k),  0,         F2(k);
       -F3(k), -F2(k), 0        ];
  sksym = dQ*Q' + Q*F*Q';
  sksymF = .5 * (sksym - sksym');   % is this OK?
  skErrormat = sksym - sksymF;
  skErrorSym(k) = norm(triu(skErrormat,1),'fro');
  skErrorDiag(k) = norm(diag(skErrormat),'fro');
  skError(k) = norm(skErrormat,'fro');
%   display(['(' mfilename ') error for skew-symmetry: ' num2str(skError(k))])
  if skError(k) > 1
    display(['skew-symmetry error: ' num2str(chg) '/' num2str(skError(k))])
%     display(sksym - sksymF)
  end
  Gammam1(k) = sksymF(1,2);
  Gammam2(k) = sksymF(2,3);
  Gammam3(k) = sksymF(1,3);
  Lambda11(k) = Lambda11(k) + sksym(1,1);
  Lambda21(k) = Lambda21(k) + sksym(1,1);
  Lambda31(k) = Lambda31(k) + sksym(1,1);
end
    
%% Save variables

output_dir = '../output/';
display(['(' mfilename ') plotting / saving outputs'])

save([output_dir 'Gammal_' sigma_name '_' num2str(resolution) '.mat'], ...
              'Gammam1','Gammam2','Gammam3','skError','skErrorSym','skErrorDiag');
% modify Lambda here!
save([output_dir 'LambdalFl_' sigma_name '_' num2str(resolution) '.mat'], ...
              'Lambda11','Lambda21','Lambda31','-append');

end
