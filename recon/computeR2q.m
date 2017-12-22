function [a,b,c,d] = computeR2q(R)
t = trace(R);
r = sqrt(1 + t);
a = .5 * r;
b = sign(R(3,2) - R(2,3))*abs(.5 * sqrt(1 + R(1,1) - R(2,2) - R(3,3)));
c = sign(R(1,3) - R(3,1))*abs(.5 * sqrt(1 - R(1,1) + R(2,2) - R(3,3)));
d = sign(R(2,1) - R(1,2))*abs(.5 * sqrt(1 - R(1,1) - R(2,2) + R(3,3)));
end