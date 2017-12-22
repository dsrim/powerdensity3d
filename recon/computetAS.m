function tAS = computetAS(Gamma,u1,u2,u3,h)

[du1dx, du1dy, du1dz] = computeGrad(u1,h);
[du2dx, du2dy, du2dz] = computeGrad(u2,h);
[du3dx, du3dy, du3dz] = computeGrad(u3,h);

N = size(u1,1);
M = numel(u1);
tAS = initMatx(3,3,N);

for j = 1:M
    S = [du1dx(j), du1dy(j), du1dz(j);
         du2dx(j), du2dy(j), du2dz(j);
         du3dx(j), du3dy(j), du3dz(j)]';
     
    gammax = getMatrixAtj(Gamma,j);
    detgamma = det(gammax);
    
    tASx = gammax/(detgamma)^(1/3) * S;
    
    tAS = setMatrixAtj(tAS,tASx,j);
    
end


end