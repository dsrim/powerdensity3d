function Tqv = Tq(q,v)
%TQ compute the map T_q v 
% q is a unit 4-vector containing coefficients in quaternions {1,e1,e2,e3}
% v is a 3-vector

% TODO check if v is a unit-vector?

v = reshape(v,3,1);     % make sure the vector is in column shape

a = q(1);   b = q(2);   c = q(3);   d = q(4);

Tq = [a^2 + b^2 - c^2 - d^2,           2*(b*c - a*d),          2*(b*d + a*c);
              2*(b*c + a*d),   a^2 - b^2 + c^2 - d^2,          2*(c*d - a*b);
              2*(b*d - a*c),           2*(c*d + a*b),  a^2 - b^2 - c^2 + d^2];

Tqv = Tq*v;

end
