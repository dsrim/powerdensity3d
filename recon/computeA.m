function A = computeA(Vijl,q,k)
%COMPUTEA
%
% computes a(q) given in eq:a at spatial location x, given by index k
%


V = Vijl(:,:,:,k);  % V evaluated at current spatial point

A = zeros(4,4,3);
a = zeros(3,3);
I = eye(3);

qbar = q;
qbar(2:4) = -qbar(2:4);

% symmetric part
Va = zeros(3,3,3);         % anti symmetric pt of V
for i = 1:3
    for j = 1:3
        Va(i,j,:) = .5*(V(i,j,:) - V(j,i,:));
    end
end

for m = 1:3
    
    em = I(:,m);
    t = Tq(qbar,em);
    
    TqbVt = zeros(3,3,3,3);    % t * Tqbar of symmetric pt of V
    for i = 1:3
    for j = 1:3
        Vsij = .5*(V(i,j,:) + V(j,i,:));
        for k = 1:3
            TqbVt(i,j,:,k) = Tq(qbar, Vsij)*t(k);
        end
    end
    end

    a(1,m) = em'*reshape(Va(2,3,:),3,1) ...
             + trace(reshape(TqbVt(3,:,2,:) - TqbVt(2,:,3,:),3,3)) ...
       + 2./3.*trace(reshape(TqbVt(2,:,:,3) - TqbVt(3,:,:,2),3,3));

    a(2,m) = em'*reshape(Va(3,1,:),3,1) ...
             + trace(reshape(TqbVt(1,:,3,:) - TqbVt(3,:,1,:),3,3)) ...
       + 2./3.*trace(reshape(TqbVt(3,:,:,1) - TqbVt(1,:,:,3),3,3));

    a(3,m) = em'*reshape(Va(1,2,:),3,1) ...
             + trace(reshape(TqbVt(2,:,1,:) - TqbVt(1,:,2,:),3,3)) ...
       + 2./3.*trace(reshape(TqbVt(1,:,:,2) - TqbVt(2,:,:,1),3,3));

    
end

% output result right multiplication by a in matrix form
for m = 1:3
    A(:,:,m) = [      0, -a(1,m), -a(2,m), -a(3,m);
                 a(1,m),       0,  a(3,m), -a(2,m);
                 a(2,m), -a(3,m),       0,  a(1,m);
                 a(3,m),  a(2,m), -a(1,m),       0];
end
end
