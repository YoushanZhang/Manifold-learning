
function out=projectTangent( base, v,method)
if nargin<3
    method = 'Sphere';
end
%% Sphere
if(strcmpi(method,'Sphere'))
    
out = v - base'* v * base;
 
 %% Kendall space
elseif (strcmpi(method,'Kendall'))
   
    
  p=base;
  num = size(p,2);
%   meanOfV = zeros(1, 2);
  meanOfV = mean(v, 2);

  hv = v;
%    Remove translation component
  for i = 1: num
    hv(:, i) = hv(:,i) - meanOfV;
  end
  
%    Project onto tangent space of sphere
  hv = hv - p * innerProduct(p, p, hv,method);
%    Horizontal projection
  skew = zeros(2, 2);
  skew(1, 2) = -1;
  skew(2, 1) = 1;

  vert = skew * p;
  out = hv - vert * innerProduct(p, hv, vert,method);

end
end
