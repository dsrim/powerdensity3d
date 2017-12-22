function Gamma = getGamma(N,gammatype)

M = N^3;


switch lower(gammatype)
  case 'id'
    
    % identity matrix
    
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);

  
  case 'scalar_smoothtori3'
    % identity matrix

    sigma_name = gammatype(8:end);
    
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);

    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    sigval = reshape(my3Dmap(sigma_name,[X(:),Y(:),Z(:)]')',1,1,M);
    Gamma = Gamma .* repmat(sigval,3,3,1);
 
  case 'aniso1'
    % gaussians located on the diagonal, sigma = .2
    
    sigma = .2;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma = zeros(3,3,M);
    Gamma(1,1,:) = 1 + gaussian(X,Y,Z,[0.5,0.0,0.0],sigma);
    Gamma(2,2,:) = 1 + gaussian(X,Y,Z,[0.0,0.5,0.0],sigma);
    Gamma(3,3,:) = 1 + gaussian(X,Y,Z,[0.0,0.0,0.5],sigma);   
    
  case 'aniso2'
    % gaussians located on the diagonal, sigma = 2
    
    sigma = 2.;
     x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma = zeros(3,3,M);
    Gamma(1,1,:) = 1 + gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(2,2,:) = 1 + gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(3,3,:) = 1 + gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
     
  case 'aniso3'
    % gaussians located on the diagonal, sigma = .5
      
    sigma = .5;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma = zeros(3,3,M);
    Gamma(1,1,:) = 1 + gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(2,2,:) = 1 + gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(3,3,:) = 1 + gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
    
  case 'aniso11only'
    % gaussian located on (1,1)-entry
      
    sigma = .5;
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma = zeros(3,3,M);
    Gamma(1,1,:) = 1 + 5*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    
  case 'aniso_xy'
    % gaussians located on (1-2,1-2) entries
    
    sigma = .5;    
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    Gamma = zeros(3,3,M);
    Gamma(1,1,:) = 2 + 2*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(1,2,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(2,1,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(2,2,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(3,3,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
      
  case 'aniso_xyz'
    % gaussians located on all entries
    
    sigma = .5;    
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
      
    Gamma(1,1,:) = 2 + 2*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(1,2,:) = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(1,3,:) = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(2,1,:) = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(2,2,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(2,3,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,1,:) = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(3,2,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,3,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);

  case 'aniso_xyz_10'
    % gaussians located on all entries
    
    sigma = 2.;    
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
      
    Gamma(1,1,:) = 2 + 10*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(1,2,:) = 2.5*gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(1,3,:) = 2.5*gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(2,1,:) = 2.5*gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(2,2,:) = 2 + 10*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(2,3,:) = 5*gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,1,:) = 5*gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(3,2,:) = 5*gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,3,:) = 2 + 10*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
      
  case 'aniso_xyz_2'
    % gaussians located on all entries
    
    sigma = 1.;    
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
      
    Gamma(1,1,:) = 2 + 2*gaussian(X,Y,Z,[0.5;0.0;0.0],sigma);
    Gamma(1,2,:) = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(1,3,:) = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(2,1,:) = gaussian(X,Y,Z,[0.0;0.0;-0.5],sigma);
    Gamma(2,2,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.5;0.0],sigma);
    Gamma(2,3,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,1,:) = gaussian(X,Y,Z,[0.0;-0.5;0.0],sigma);
    Gamma(3,2,:) = gaussian(X,Y,Z,[-0.5;0.0;0.0],sigma);
    Gamma(3,3,:) = 2 + 2*gaussian(X,Y,Z,[0.0;0.0;0.5],sigma);
 

  case 'aniso_phi'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    M = N^3;
    Gamma = zeros(3,3,M);
    
    phi = [ 1; -1; 0]/sqrt(2);
    matphi = phi*phi';
    lam1 = lam(X,Y,Z);
    
    Gamma(1,1,:) = lam1*matphi(1,1) + 1.;
    Gamma(1,2,:) = lam1*matphi(1,2);
    Gamma(1,3,:) = lam1*matphi(1,3);
    Gamma(2,1,:) = lam1*matphi(2,1);
    Gamma(2,2,:) = lam1*matphi(2,2) + 1.;
    Gamma(2,3,:) = lam1*matphi(2,3);
    Gamma(3,1,:) = lam1*matphi(3,1);
    Gamma(3,2,:) = lam1*matphi(3,2);
    Gamma(3,3,:) = lam1*matphi(3,3) + 1.;
   
    
   case 'aniso_tori'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    k = 1;              % max. conductivity for tori
    r = 0.1;            % small radius of tori
    R = 0.8;            % large radius of tori
       
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
       
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2

    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
  case 'aniso_tori_2'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    k = 2;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                  
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2

    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

  case 'aniso_tori_4'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    k = 4;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
     
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);

    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
  case 'aniso_tori_6'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
    
    k = 6;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

  case 'aniso_tori_8'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 8;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

  case 'aniso_tori_10'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 10;             % max. conductivity for tori
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
  
    k = 10;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

  case 'aniso_tori_11'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 11;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);


  case 'aniso_tori_12'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 12;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);



  case 'aniso_tori_15'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 15;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

  case 'aniso_tori_25'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 25;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);


  case 'aniso_tori_25'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
  
    k = 25;             % max. conductivity for tori
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
  
    k = 50;             % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [0; 0; -0.5];  % center of tori 1
    phi1 = [1; 0; 0];   % direction for tori 1
                                                   
    x2 = [0; 0; 0.5];   % center of tori 2
    phi2 = [0; 1; 0];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);

      
  case 'aniso_tori2_2'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
                
    k = 2;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [-0.5; 0; 0];  % center of tori 1
    phi1 = [0; 1; 0];   % direction for tori 1
                                               
    x2 = [0.5; 0; 0];   % center of tori 2
    phi2 = [0; 0; 1];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
        
  case 'aniso_tori2_4'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);
            
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
            
    k = 6;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
        
    x1 = [-0.5; 0; 0];  % center of tori 1
    phi1 = [0; 1; 0];   % direction for tori 1
                                                   
    x2 = [0.5; 0; 0];   % center of tori 2
    phi2 = [0; 0; 1];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
  case 'aniso_tori2_8'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);

    k = 8;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
    
    x1 = [-0.5; 0; 0];  % center of tori 1
    phi1 = [0; 1; 0];   % direction for tori 1
                                               
    x2 = [0.5; 0; 0];   % center of tori 2
    phi2 = [0; 0; 1];   % direction for tori 2
        
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
        
        
  case 'aniso_tori2_10'
        
    x = linspace(-1,1,N);
    [X,Y,Z] = ndgrid(x,x,x);

    k = 10;              % max. conductivity for tori
    r = 0.1;            % small radius of tori          
    R = 0.8;            % large radius of tori
    
    x1 = [-0.5; 0; 0];  % center of tori 1
    phi1 = [0; 1; 0];   % direction for tori 1
                                               
    x2 = [0.5; 0; 0];   % center of tori 2
    phi2 = [0; 0; 1];   % direction for tori 2
    
    Gamma = zeros(3,3,M);    
    Gamma(1,1,:) = ones(M,1);
    Gamma(2,2,:) = ones(M,1);
    Gamma(3,3,:) = ones(M,1);
    Gamma = Gamma + anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R);
end

end


function Gamma = anisoTori(N,X,Y,Z,x1,phi1,x2,phi2,k,r,R)

M = N^3;
Gamma = zeros(3,3,M);
        
X0 = [X(:)';Y(:)';Z(:)'];
        
X1 = X0 - repmat(x1,1,M);
X2 = X0 - repmat(x2,1,M);

   
P1 = X1 - phi1*(phi1'*X1);
P1 = R*P1 ./ repmat(norm3(P1),3,1) + repmat(x1,1,M);
Chi1 = reshape(exp( -.5*norm3(X0 - P1).^2 / r^2),1,1,M);
Y1 = cross(X1./ repmat(norm3(X1),3,1),repmat(phi1,1,M),1);

Gamma1 =  repmat(Chi1,3,3,1) .* outer3(Y1,Y1);


P2 = X2 - phi2*(phi2'*X2);
P2 = R*P2 ./ repmat(norm3(P2),3,1) + repmat(x2,1,M);
Chi2 = reshape(exp( -.5*norm3(X0 - P2).^2 / r^2),1,1,M);
Y2 = cross(X2./ repmat(norm3(X2),3,1),repmat(phi2,1,M),1);

Gamma2 =  repmat(Chi2,3,3,1) .* outer3(Y2,Y2);

Gamma = k*Gamma1 + k*Gamma2;


end

function val = gaussian(X,Y,Z,X0,sigma)
% the output will be a vector of size (N^3,1)

val = exp(-((X(:)-X0(1)).^2 +(Y(:)-X0(2)).^2 +(Z(:)-X0(3)).^2)/(2*sigma^2));

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

val =  5. * gaussian(X,Y,Z,[.5,.5,0.],.1);


end


