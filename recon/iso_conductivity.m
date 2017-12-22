function [fval,f] = iso_conductivity(type,const,coords)
%ISO_CONDUCTIVITY
% iso_conductivity(type,coords)
% 
% type:    string designating conductivity type
% coords:  an 3 by N array, e.g., [X(:),Y(:),Z(:)]'
% const:   parameter >=0 amount of perturbation from one
% 
% evaluates the conductivity at the grid points
% then returns a function handle g: for use with PDE toolbox
%
% prepending 'log' to the type string will return log(conductivity).
% of the given type.
%
% supported type:
%       'balls'         characteristic function of scattered balls
%       'tinyballs'     characteristic function of smaller balls
%       'gaussians'     scattered smooth gaussians function
%       'tori'          charactersitic function of tori
%       'smoothtori'    smooth gaussian tori


if strcmp(type(1:3),'log')
  logswitch = true;
  type = type(4:end);
else
  logswitch = false;
end


switch type
    case 'balls'
        
        s = 0.2;    % width of the ball
        f = @(coords) elem('one', coords) ...
             + const*(elem('ball', coords, s,  0.5,  0.5,  0.5) ...
                    + elem('ball', coords, s,  0.5,  0.5, -0.5) ...
                    + elem('ball', coords, s,  0.5, -0.5,  0.5) ...
                    + elem('ball', coords, s, -0.5, -0.5, -0.5));
        
    case 'tinyballs'
        
        s = 0.05;   % width of the ball
        f = @(coords) elem('one', coords) ...
             + const*(elem('ball', coords, s,  0.5,  0.5,  0.5) ...
                    + elem('ball', coords, s,  0.5,  0.5, -0.5) ...
                    + elem('ball', coords, s,  0.5, -0.5,  0.5) ...
                    + elem('ball', coords, s, -0.5, -0.5, -0.5) );
        
    case 'gaussians'

        s = 0.2;    % width of the gaussians
        f = @(coords) elem('one', coords) ...
             + const*(elem('gaussian', coords, s,  0.5,  0.5,  0.5) ...
                    + elem('gaussian', coords, s,  0.5,  0.5, -0.5) ...
                    + elem('gaussian', coords, s,  0.5, -0.5,  0.5) ...
                    + elem('gaussian', coords, s, -0.5, -0.5, -0.5) ...
                    + elem('gaussian', coords, s, -0.5,  0.5,  0.5) );
                     
    case 'tori'
        
        R = 0.4;    % radius of the torus
        s = 0.05;   % width of the torus
        f = @(coords) elem('one', coords) ...
            + const*(double(torus(R,0,0,-0.2,1,0,0,coords) <= s) ...
                   + double(torus(R,0,0, 0.2,0,1,0,coords) <= s) );

    case 'smoothtori'    
        
        R = 0.8;    % radius of the torus
        s = 0.05;   % width of the torus
        f = @(coords) elem('one', coords) ...
            + const*(exp(-(torus(R,0,0,-0.5,1,0,0,coords)).^2/(2*s^2)) ...
                   + exp(-(torus(R,0,0, 0.5,0,1,0,coords)).^2/(2*s^2)) );  
        
    otherwise 
        
        f = @(coords) elem('one', coords);
        
end

fval = f(coords);

if logswitch
    fval = log(fval);
    
    f = @(x) log(f(x));
end
end



function fval = elem(type,coords,sig,x,y,z)

switch type
    case 'one'
        fval = ones(1,size(coords,2));
    case 'ball'
        fval = double( (coords(1,:)-x).^2 + (coords(2,:)-y).^2 ...
            + (coords(3,:)-z).^2 <= sig^2);
    case 'gaussian'
        fval = exp(-((coords(1,:)-x).^2 + (coords(2,:)-y).^2 ...
            + (coords(3,:)-z).^2)/(2*sig^2));
    otherwise 
        fval = zeros(size(coords,1));
end

end

function fval = torus(R,x0,y0,z0,nx,ny,nz,coords)
% distance function to a circle of center (x0,y0,z0), radius R and normal 
% nx,ny,nz

n = norm([nx,ny,nz],2); 
nx = nx/n; ny = ny/n; nz = nz/n;    % normalize vector [nx,ny,nz] 
fval = zeros(1,size(coords,2));

dx = coords(1,:)-x0;                
dy = coords(2,:)-y0;
dz = coords(3,:)-z0;

dXdotN = nx*dx + ny*dy + nz*dz;

dirx = dx - dXdotN*nx;
diry = dy - dXdotN*ny;
dirz = dz - dXdotN*nz;

normdir = dirx.*dirx + diry.*diry + dirz.*dirz;
idx = (abs(normdir) < 1e-4);

fval(idx) = sqrt(R^2 + dx(idx).^2 + dy(idx).^2 + dz(idx).^2);

Cx = x0 + R*dirx(~idx)./sqrt(normdir(~idx));
Cy = y0 + R*diry(~idx)./sqrt(normdir(~idx));
Cz = z0 + R*dirz(~idx)./sqrt(normdir(~idx));

fval(~idx) = sqrt( (Cx-coords(1,~idx)).^2 + ...
                   (Cy-coords(2,~idx)).^2 + ...
                   (Cz-coords(3,~idx)).^2 );

end

    
