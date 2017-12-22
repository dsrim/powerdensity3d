function Gamma = getGamma_struct(N,gammatype)


switch lower(gammatype)
  case 'id'
    
    Gamma.m11 = ones(N,N,N);
    Gamma.m12 = zeros(N,N,N);
    Gamma.m13 = zeros(N,N,N);
    Gamma.m21 = zeros(N,N,N);
    Gamma.m22 = ones(N,N,N);
    Gamma.m23 = zeros(N,N,N);
    Gamma.m31 = zeros(N,N,N);
    Gamma.m32 = zeros(N,N,N);
    Gamma.m33 = ones(N,N,N);
    
  case 'aniso1'
    
    sigma = .2;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma.m11 = 1 + gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma.m12 = zeros(N,N,N);
    Gamma.m13 = zeros(N,N,N);
    Gamma.m21 = zeros(N,N,N);
    Gamma.m22 = 1 + gaussian(X,Y,Z,[0.;0.5;0.0],sigma);
    Gamma.m23 = zeros(N,N,N);
    Gamma.m31 = zeros(N,N,N);
    Gamma.m32 = zeros(N,N,N);
    Gamma.m33 = 1 + gaussian(X,Y,Z,[0.;0.0;0.5],sigma);   
    
  case 'aniso2'
      
    sigma = 2.;
     x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma.m11 = 1 + gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma.m12 = zeros(N,N,N);
    Gamma.m13 = zeros(N,N,N);
    Gamma.m21 = zeros(N,N,N);
    Gamma.m22 = 1 + gaussian(X,Y,Z,[0.;0.5;0.0],sigma);
    Gamma.m23 = zeros(N,N,N);
    Gamma.m31 = zeros(N,N,N);
    Gamma.m32 = zeros(N,N,N);
    Gamma.m33 = 1 + gaussian(X,Y,Z,[0.;0.0;0.5],sigma);
     
    
  case 'aniso3'
      
    sigma = .5;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma.m11 = 1 + gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma.m12 = zeros(N,N,N);
    Gamma.m13 = zeros(N,N,N);
    Gamma.m21 = zeros(N,N,N);
    Gamma.m22 = 1 + gaussian(X,Y,Z,[0.;0.5;0.0],sigma);
    Gamma.m23 = zeros(N,N,N);
    Gamma.m31 = zeros(N,N,N);
    Gamma.m32 = zeros(N,N,N);
    Gamma.m33 = 1 + gaussian(X,Y,Z,[0.;0.0;0.5],sigma);
    
    
  case 'aniso11only'
      
    sigma = .5;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma.m11 = 1 + 5*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma.m12 = zeros(N,N,N);
    Gamma.m13 = zeros(N,N,N);
    Gamma.m21 = zeros(N,N,N);
    Gamma.m22 = ones(N,N,N);
    Gamma.m23 = zeros(N,N,N);
    Gamma.m31 = zeros(N,N,N);
    Gamma.m32 = zeros(N,N,N);
    Gamma.m33 = ones(N,N,N);
    
    case 'aniso_xy'
        sigma = .5;    
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        
        Gamma.m11 = 2 + 2*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
        Gamma.m12 = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
        Gamma.m13 = zeros(N,N,N);
        Gamma.m21 = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
        Gamma.m22 = 2 + 2*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
        Gamma.m23 = zeros(N,N,N);
        Gamma.m31 = zeros(N,N,N);
        Gamma.m32 = zeros(N,N,N);
        Gamma.m33 = 2 + 2*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
    
    case 'aniso_xyz'
        sigma = .5;    
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        
        Gamma.m11 = 2 + 2*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
        Gamma.m12 = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
        Gamma.m13 = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
        Gamma.m21 = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
        Gamma.m22 = 2 + 2*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
        Gamma.m23 = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
        Gamma.m31 = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
        Gamma.m32 = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
        Gamma.m33 = 2 + 2*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
     
    case 'aniso_phi'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        Gamma = initMatx(3,3,N);
        
        for j = 1:M
            gammax = eye(3,3) ...
                    + lam(X(j),Y(j),Z(j))*phi(X,Y,Z)*phi(X,Y,Z)';
            Gamma = setMatrixAtj(Gamma,gammax,j);
        end

   case 'aniso_tori'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 1;              % max. conductivity for tori
        r = 0.1;            % small radius of tori
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
        
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2

        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
    case 'aniso_tori_2'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 2;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

    case 'aniso_tori_4'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 4;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
    case 'aniso_tori_6'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 6;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

    case 'aniso_tori_8'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 8;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

    case 'aniso_tori_10'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 10;             % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

    case 'aniso_tori_50'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 50;             % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of tori 1
        phi1 = [1; 0; 0];   % direction for tori 1
                                                   
        x2 = [0; 0; 0.5];   % center of tori 2
        phi2 = [0; 1; 0];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);


        
    case 'aniso_tori2_2'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 2;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [-0.5; 0; 0];  % center of tori 1
        phi1 = [0; 1; 0];   % direction for tori 1
                                                   
        x2 = [0.5; 0; 0];   % center of tori 2
        phi2 = [0; 0; 1];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
        
    case 'aniso_tori2_4'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 4;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [-0.5; 0; 0];  % center of tori 1
        phi1 = [0; 1; 0];   % direction for tori 1
                                                   
        x2 = [0.5; 0; 0];   % center of tori 2
        phi2 = [0; 0; 1];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
        
    case 'aniso_tori2_6'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 6;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [-0.5; 0; 0];  % center of tori 1
        phi1 = [0; 1; 0];   % direction for tori 1
                                                   
        x2 = [0.5; 0; 0];   % center of tori 2
        phi2 = [0; 0; 1];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
    case 'aniso_tori2_8'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 8;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [-0.5; 0; 0];  % center of tori 1
        phi1 = [0; 1; 0];   % direction for tori 1
                                                   
        x2 = [0.5; 0; 0];   % center of tori 2
        phi2 = [0; 0; 1];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
        
    case 'aniso_tori2_10'
        
        x = linspace(-1,1,N);
        [X,Y,Z] = ndgrid(x,x,x);
        M = N^3;
        k = 10;              % max. conductivity for tori
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [-0.5; 0; 0];  % center of tori 1
        phi1 = [0; 1; 0];   % direction for tori 1
                                                   
        x2 = [0.5; 0; 0];   % center of tori 2
        phi2 = [0; 0; 1];   % direction for tori 2
        
        Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
end

end


function Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R)

        Gamma = initMatx(3,3,N);
        M = N^3;

        for j = 1:M
            xx = [X(j); Y(j); Z(j)];
            
            P1 = (xx - x1) - (xx - x1)'*phi1*phi1;
            P1 = R*P1/norm(P1) + x1;
            chi1 = exp( -norm(xx- P1,2)^2/(2*r^2));
            yy1 = cross((xx - x1)/norm(xx - x1),phi1);
            gamma1 = chi1*(yy1*yy1');
%             gamma1 = chi1*eye(3);
            
            P2 = (eye(3) - phi2*phi2')*(xx - x2);
            P2 = R*P2/norm(P2) + x2;
            
            chi2 = exp( -norm(xx- P2,2)^2/ (2*r^2));
            yy2 = cross((xx - x2)/norm(xx - x2),phi2);
            gamma2 = chi2*(yy2*yy2');
%             gamma2 = chi2*eye(3);
            
            gammax = eye(3) + k*gamma1 + k*gamma2;
            
            Gamma = setMatrixAtj(Gamma,gammax,j);
        end       

end

function val = gaussian(X,Y,Z,X0,sigma)
% sigma = .2;

N = size(X,1);
M = numel(X);
val = zeros(N,N,N);

for j = 1:M
   xx = [X(j);Y(j);Z(j)];
   val(j) = exp(- norm(xx - X0)^2/(2*sigma^2));
end

end



function val = phi(X,Y,Z)

% al = 0;
% be = pi/4;
% 
% rotA = [cos(al) -sin(al) 0 ;
%         sin(al)  cos(al) 0 ;
%            0  0  1];
% v  = [cos(be);
%       sin(be);
%       0];
% 
% val = rotA*v;


val = [ 0; 1; 1]/sqrt(2);


end

function val = lam(X,Y,Z)

val =  5. * gaussian(X,Y,Z,[0.;.5;0.],.2);


end


