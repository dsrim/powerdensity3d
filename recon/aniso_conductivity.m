function c = aniso_conductivity(type,const)

switch type
    
    case 'anisophi' % 1D anisotropic perturbation
        phi = [1; 1; 0]/sqrt(2);
        
        % scalar blob centered at zero
        sig2 = (0.1)^2;
        cg = @(x,y,z) exp(-(x.^2 + y.^2 + z.^2)/(2*sig2));
        
        c = @(x,y,z) reshape(eye(3),9,1) * ones(size(x)) ...
            + const .* reshape(phi*phi',9,1) * cg(x,y,z); 
        
    case 'isotori'
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of torus 1
        phi1 = [1; 0; 0];   % direction for torus 1
                                                   
        x2 = [0; 0; 0.5];   % center of torus 2
        phi2 = [0; 1; 0];   % direction for torus 2
        
        P1 = @(x,y,z) (eye(3)-phi1*phi1')*([x;y;z]-x1*ones(size(x)));
        Q1 = @(x,y,z) R * P1(x,y,z)./([1;1;1]*sqrt(sum(P1(x,y,z).^2,1))) + x1*ones(size(x));
        chi1 = @(x,y,z) exp( - sum(([x;y;z] - Q1(x,y,z)).^2, 1)/(2*r^2));
        
        P2 = @(x,y,z) (eye(3)-phi2*phi2')*([x;y;z]-x2*ones(size(x)));
        Q2 = @(x,y,z) R * P2(x,y,z)./([1;1;1]*sqrt(sum(P2(x,y,z).^2,1))) + x2*ones(size(x));
        chi2 = @(x,y,z) exp( - sum(([x;y;z] - Q2(x,y,z)).^2, 1)/(2*r^2));
        
        c = @(x,y,z) reshape(eye(3),9,1) * (ones(size(x)) + const * (chi1(x,y,z) + chi2(x,y,z)));            
  
        
    case 'anisotori'
        r = 0.1;            % small radius of tori          
        R = 0.8;            % large radius of tori
        
        x1 = [0; 0; -0.5];  % center of torus 1
        phi1 = [1; 0; 0];   % direction for torus 1
                                                   
        x2 = [0; 0; 0.5];   % center of torus 2
        phi2 = [0; 1; 0];   % direction for torus 2
        
        P1 = @(x,y,z) (eye(3)-phi1*phi1')*([x;y;z]-x1*ones(size(x)));
        Q1 = @(x,y,z) R * P1(x,y,z)./([1;1;1]*sqrt(sum(P1(x,y,z).^2,1))) + x1*ones(size(x));
        n1 = @(x,y,z) sqrt(sum( ([x;y;z]-x1*ones(size(x))).^2, 1));
        chi1 = @(x,y,z) exp( - sum(([x;y;z] - Q1(x,y,z)).^2, 1)/(2*r^2));
        yy1 = @(x,y,z) [phi1(3)*(y-x1(2)) - phi1(2)*(z-x1(3)); 
                        phi1(1)*(z-x1(3)) - phi1(3)*(x-x1(1));
                        phi1(2)*(x-x1(1)) - phi1(1)*(y-x1(2))] ./ ([1;1;1]*n1(x,y,z));
        
        
        P2 = @(x,y,z) (eye(3)-phi2*phi2')*([x;y;z]-x2*ones(size(x)));
        Q2 = @(x,y,z) R * P2(x,y,z)./([1;1;1]*sqrt(sum(P2(x,y,z).^2,1))) + x2*ones(size(x));
        n2 = @(x,y,z) sqrt(sum( ([x;y;z]-x2*ones(size(x))).^2, 1));
        chi2 = @(x,y,z) exp( - sum(([x;y;z] - Q2(x,y,z)).^2, 1)/(2*r^2));
        yy2 = @(x,y,z) [phi2(3)*(y-x2(2)) - phi2(2)*(z-x2(3)); 
                        phi2(1)*(z-x2(3)) - phi2(3)*(x-x2(1));
                        phi2(2)*(x-x2(1)) - phi2(1)*(y-x2(2))] ./ ([1;1;1]*n2(x,y,z));
        
        
        c1 = @(x,y,z) (ones(9,1)*chi1(x,y,z)) .* kron([1;1;1],yy1(x,y,z)) .* kron(yy1(x,y,z),[1;1;1]);
        c2 = @(x,y,z) (ones(9,1)*chi2(x,y,z)) .* kron([1;1;1],yy2(x,y,z)) .* kron(yy2(x,y,z),[1;1;1]);
                
        c = @(x,y,z) reshape(eye(3),9,1)*ones(size(x)) ...
            + const * (c1(x,y,z) + c2(x,y,z));        
    
    otherwise % constant identity matrix
        c = @(x,y,z) reshape(eye(3),9,1)*ones(size(x)); 
end


end
