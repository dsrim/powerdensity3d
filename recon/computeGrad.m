function [dudx, dudy, dudz] = computeGrad(u,h,varargin)
% computes the gradient of u
%
% computeGrad(u,h)
%
% u is a (N,N,N) matrix on uniform grid with meshgrid ordering
% h is the width of the uniform grid
%

s = size(u);
if (length(s) == 3) && (s(1) == s(2)) && (s(2) == s(3))
  returnvector = false;  
else
  M = length(u);
  N = floor(M^(1/3));
  if N^3 ~= M
    N = N+1;
  end
  u = reshape(u,N,N,N);
  
  returnvector = true;
end


if ~isempty(varargin) && strcmp(varargin{1}, 'meshgrid')
    
    % use meshgrid ordering
    
    dudx = (u(:,[2:end,end],:) - u(:,[1,1:(end-1)],:))/(2*h);
    dudx(:,1,:) = (-u(:,3,:) + 4*u(:,2,:) - 3*u(:,1,:))/(2*h);
    dudx(:,end,:) = (3*u(:,end,:) - 4*u(:,end-1,:) + u(:,end-2,:))/(2*h);

    dudy = (u([2:end,end],:,:) - u([1,1:(end-1)],:,:))/(2*h);
    dudy(1,:,:) = (-u(3,:,:) + 4*u(2,:,:) - 3*u(1,:,:))/(2*h);
    dudy(end,:,:) = (3*u(end,:,:) - 4*u(end-1,:,:) + u(end-2,:,:))/(2*h);

    dudz = (u(:,:,[2:end,end]) - u(:,:,[1,1:(end-1)]))/(2*h);
    dudz(:,:,1) = (-u(:,:,3) + 4*u(:,:,2) - 3*u(:,:,1))/(2*h);
    dudz(:,:,end) = (3*u(:,:,end) - 4*u(:,:,end-1) + u(:,:,end-2))/(2*h);
else
    
    % use ndgrid ordering
    
    dudx = (u([2:end,end],:,:) - u([1,1:(end-1)],:,:))/(2*h);
    dudx(1,:,:) = (-u(3,:,:) + 4*u(2,:,:) - 3*u(1,:,:))/(2*h);
    dudx(end,:,:) = (3*u(end,:,:) - 4*u(end-1,:,:) + u(end-2,:,:))/(2*h);
    
    dudy = (u(:,[2:end,end],:) - u(:,[1,1:(end-1)],:))/(2*h);
    dudy(:,1,:) = (-u(:,3,:) + 4*u(:,2,:) - 3*u(:,1,:))/(2*h);
    dudy(:,end,:) = (3*u(:,end,:) - 4*u(:,end-1,:) + u(:,end-2,:))/(2*h);

    dudz = (u(:,:,[2:end,end]) - u(:,:,[1,1:(end-1)]))/(2*h);
    dudz(:,:,1) = (-u(:,:,3) + 4*u(:,:,2) - 3*u(:,:,1))/(2*h);
    dudz(:,:,end) = (3*u(:,:,end) - 4*u(:,:,end-1) + u(:,:,end-2))/(2*h);
    
end


if (returnvector == 1)
  % if received a vector, return a vector
  
  dudx = dudx(:);
  dudy = dudy(:);
  dudz = dudz(:);
end

end
