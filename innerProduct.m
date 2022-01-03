function result= innerProduct( base, v, w,method)
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%% 
if nargin<4
    method = 'Sphere';
end
%% Sphere
if(strcmpi(method,'Sphere'))
    
dimension=size(base);
 result = 0.0;
 for i=1:1: dimension
    result = result+v(i) * w(i);
 end
 
 
 %% Kendall space
elseif (strcmpi(method,'Kendall'))
  result=sum(dot(v,w));

  
  %%  Grassmannian manifold
elseif (strcmpi(method,'Grassmannian'))
    
    result = trace(v' * w);
end

  
end
 
