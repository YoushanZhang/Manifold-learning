function [pi, mu, W, Lambda, sigma2,clusters] = k_means_PGA(t, p, method, variance_level)
 %% This function performs K means principal geodesic analysis (PGA) on manifolds for muti clusters
% Input: 
%   t: d-by-n matrix
%   p: number of clusters
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
if nargin<3
    variance_level = 1;
    method='Sphere';
end

if size(t,3)==1    
   [d, Npoints] = size(t);
else
   [d,~, Npoints] = size(t); 
end

%% k-means

% initialization
init_centers = randi(Npoints, p, 1);
if size(t,3)>1
    mu = t(:, :, init_centers);
else
    mu = t(:, init_centers);
end
distance_square = zeros(Npoints, p);
clusters = zeros(Npoints, 1);
Dold = -2;
D = -1;

while Dold ~= D
    Dold = D;
    
    if size(t,3)>1
        for c = 1:p
            distance_square(:,c) = zeros(Npoints, 1);
         for j=1:1:Npoints
             distance_square(j,c) = distance_square(j,c) +...
                 (norm(log_Map(t(:,:,j),mu(:,:,c),method)).^2)';
         end
        end
    else        
    % assign classes
    for c = 1:p
        distance_square(:,c) = zeros(Npoints, 1);
     for j=1:1:Npoints
        distance_square(j,c) = distance_square(j,c) + (norm(log_Map(t(:,j),mu(:,c),method)).^2)';
     end
    end
    end
    [distmin, clusters] = min(distance_square, [], 2);
% distance_square
%     clusters
    % compute distortion
    D = sum(distmin);
    % compute new centers
    if size(t,3)>1
       for c=1:p
        mu(:, :,c) = intrinsicMean(t(:,:, clusters == c),method);    
       end
    else      
       for c=1:p
        mu(:, c) = intrinsicMean(t(:, clusters == c),method);    
       end
    end
end



%% parameter initialization
pi = zeros(p, 1);
if size(t,3)==1
   W = zeros(size(t,1),size(t,1),p);
   Lambda=zeros(size(t,1),size(t,1),p);
else
   W = zeros(size(t,1)*size(t,2),size(t,1)*size(t,2),p);
   Lambda=zeros(size(t,1)*size(t,2),size(t,1)*size(t,2),p);
end
sigma2 = zeros(p, 1);
for c=1:p
%     W{c} = normrnd(0, variance_level, [d, q]);
    I = find(clusters==c);
%     size(t(:,I),2)
    if size(t,3)==1
       [~,Lambda(:,:,c),W(:,:,c)] = PGA(t(:,I),method);
    else
       [~,Lambda(:,:,c),W(:,:,c)] = PGA(t(:,:,I),method);
    end
    
    pi(c) = numel(I)/Npoints;
    if (nargin <5)
        sigma2(c) = mean(distmin(I))/d;
    else
        sigma2(c)= abs(normrnd(0,variance_level/10));
    end
end

end
