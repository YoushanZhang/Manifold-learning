function out=log_Map(base,p,method)
% logMap is the inverse of the exponential map. It returns the initial
% veloctiy vector for the geodesic segment between base and p.
%Input:
%---base: a base point on a manifold
%---p: end point on a manifold
%---method: select one of three manifolds: sphere, Kendall, and Grassmannian
%Output:initial veloctiy vector
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%  
if nargin<3
    method = 'Sphere';
end


if(strcmpi(method,'Sphere'))
     cosTheta = innerProduct(base, base, p);
     t = p - base * cosTheta;
     length = norm(t);

  if(length < 1.0e-12 | cosTheta >= 1.0 | cosTheta <= -1.0)
    out= zeros(size(t));
  else
    out = t * (acos(cosTheta) / length);
  end
%   

elseif (strcmpi(method,'Kendall'))
    %% Kendall space size should be 2 * N
 
  q=p;p=base;
  m = p * q';
  [U, s, V]=svd(m,'econ');
  rotation = U * V';
  qRot = rotation * q;
  cosTheta = innerProduct(p, p, qRot, method);
  t = qRot - p * cosTheta;
  length = sqrt(normSquared(p, t, method));
  if(length < 1.0e-12 | cosTheta >= 1.0 | cosTheta <= -1.0)
%     out= zeros(size(t));
    out= p-q; % to avoid the NaN in svd
  else
  out= t * (acos(cosTheta) / length);
  end

  
elseif (strcmpi(method,'Grassmannian'))  

  dim = size(base,1);
  I=eye(dim, dim);
  temp1 = I - base * base';
  temp2 = base' * p;
  X = temp1 * p / temp2;
  [U, s, V]=svd(X,'econ');
   theta = atan(s);

  out= U * theta * V';
end
  
  
end
