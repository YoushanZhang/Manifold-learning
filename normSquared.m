function result= normSquared(p,v,method)
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%
if nargin<3
    method = 'Sphere';
end
%% Sphere
if(strcmpi(method,'Sphere'))
    
  result=innerProduct(p, v, v);
  
  %% Kendall space
elseif (strcmpi(method,'Kendall'))
  
  result=innerProduct(p, v, v,'Kendall'); 
  
%%  Grassmannian manifold
elseif (strcmpi(method,'Grassmannian'))
  
  result=innerProduct(p, v, v,'Grassmannian'); 
  
end

end