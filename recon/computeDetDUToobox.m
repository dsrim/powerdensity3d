function detDU = computeDetDUToobox(X,Y,Z,u1,u2,u3)

[u1x,u1y,u1z] = u1.evaluateGradient(X,Y,Z);
[u2x,u2y,u2z] = u2.evaluateGradient(X,Y,Z);
[u3x,u3y,u3z] = u3.evaluateGradient(X,Y,Z);

detDU = u1x.*(u2y.*u3z-u3y.*u2z) ...
       - u2x.*(u1y.*u3z-u1z.*u3y) ...
       + u3x.*(u1y.*u2z-u1z.*u2y);
   
end