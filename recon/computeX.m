function X = computeX(curve,Vijl,Ball)
%COMPUTEX
%
% solve the quaternionic system
%
%     q0
% d   q1
% dt  q2
%     q3
%

N = length(curve);       % grid size can be inferred from size of curve

X = zeros(4,N);         % initialize variable
X(:,1) = [1;0;0;0];     % set initial condition as identity
h = 2/(N-1);            % grid width
I = eye(4);

for j = 1:(N-1)
    
    k = curve(j,1);      % index for current point
    v = curve(j,2:4);    % integrating direction in (x,y,z) coordinates
    
    Ak = computeA(Vijl,X(:,j),k);
    Bk = Ball(:,:,k,:);

    A = v(1)*Ak(:,:,1) + v(2)*Ak(:,:,2) + v(3)*Ak(:,:,3);
    B = v(1)*Bk(:,:,1) + v(2)*Bk(:,:,2) + v(3)*Bk(:,:,3);

    anorm = norm(A(1,:));
    bnorm = norm(B(1,:));

    X(:,j+1) =  (cos(.5*h*anorm)*I + sin(.5*h*anorm)*A/anorm) ...
               *(cos(.5*h*bnorm)*I + sin(.5*h*bnorm)*B/bnorm) * X(:,j);
    
end

end
