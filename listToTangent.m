%% manifolds
function t=listToTangent(l,method)

if nargin<4
    method = 'Sphere';
end
%% Sphere
if(strcmpi(method,'Sphere'))
    
t=l;

 %% Kendall space
elseif (strcmpi(method,'Kendall'))

t=reshape(l,[2, length(l)/2]);
   
end


