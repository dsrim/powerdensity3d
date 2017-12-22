function val = computeDetU(u1,u2,u3)


N = size(u1,1);
M = N^3;
h = 2./(N-1);

val = zeros(N,N,N);

[du1dx, du1dy, du1dz] = computeGrad(u1,h);
[du2dx, du2dy, du2dz] = computeGrad(u2,h);
[du3dx, du3dy, du3dz] = computeGrad(u3,h);

for j=1:M
    A = [du1dx(j), du1dy(j), du1dz(j);
         du2dx(j), du2dy(j), du2dz(j);
         du3dx(j), du3dy(j), du3dz(j)];

    val(j) = det(A);

end

