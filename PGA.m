 function [mu, lambda, vector]=PGA(data,method)
 %% This function performs principal geodesic analysis (PGA) on manifolds
% Input: 
%   data: d-by-n matrix
%   method: Sphere, Kendall, Grassmannian
% Output:
%   mu: mean center
%   lambda: EigenValues
%   vector: Eigenvectors
%% Created by Youshan Zhang
%  First version: March 2018
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%  
if nargin<2
    method = 'Sphere';
end
 
if size(data,3)>1
   c=size(data,3);
else
   c=size(data,2);
end

 mu=intrinsicMean(data,method);
sum=0;
if size(data,3)>1
    for i=1:1:c
       mu_each=log_Map(mu, data(:,:,i),method);
       sum=sum+reshape(mu_each,[],1)*reshape(mu_each,1,[]);
    end
else
    for i=1:1:c  
       mu_each=log_Map(mu, data(:,i),method); 
       sum=sum+mu_each*mu_each';
    end
end

S= sum/c;

[vector,lambda, ~]=svd(S);


 end
