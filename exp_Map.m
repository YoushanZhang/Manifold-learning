function out = exp_Map(base,t,method)
% manifold exponential map
% expMap gives the Riemannian exponential map (i.e., geodesic segment
% starting at the base with initial velocity v).
%Input:
%---base: a base point on a manifold
%---t: time t
%---method: select one of three manifolds: sphere, Kendall, and Grassmannian
%Output:
%---out: Another point on manifold give base and t
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%
if nargin<3
    method = 'Sphere';
end

if(strcmpi(method,'Sphere'))
%   theta=norm(t);
  theta = sqrt(normSquared(base, t));
  if(theta < 1.0e-12)
    out=base;
  else
      
  p = cos(theta) * base + (sin(theta) / theta) * t;
  pNorm = sqrt(normSquared(p, p));
%   out= p / norm(p);
  out= p / pNorm;
  end
elseif (strcmpi(method,'Kendall')) 
    
 theta = sqrt(normSquared(base, t,method));
     
 if(theta < 1.0e-12)
    out= base;
 else
      
  q = cos(theta) * base + (sin(theta) / theta) * t;
  out= q / sqrt(normSquared(q, q,method));
  
 end
elseif (strcmpi(method,'Grassmannian'))  
    
   normt =  sqrt(normSquared(t, t,method));

  if(normt < 1.0e-12)
    out= base;
  end

  %  Compact SVD
  [U, s, V]=svd(t,'econ');

  target = (base * V  * diag(cos(diag(s))) + U * diag(sin(diag(s)))) * V';

  out= target;
    
end
end



