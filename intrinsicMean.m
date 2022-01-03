function m = intrinsicMean(xi,method)
%
% return the intrinsic mean of xi
%
if nargin<2
    method = 'Sphere';
end
    
step = 1.0; % can be changed during iteration
maxIter = 400;

if size(xi,3)==1
   N = size(xi,2); % number of points
   m = xi(:,1); % initial guess
else
   N = size(xi,3); % number of points
   m = xi(:,:,1); % initial guess
end
i = 1;
tol= 1.0000e-12;
while i <= maxIter
    pm = m;
    
    % iterate
    ss = zeros(size(xi));
    if size(xi,3)>1
        for j = 1:N   
         ss(:,:,j) = log_Map(m,xi(:,:,j),method);   
        end
    else
        for j = 1:N
         ss(:,j) = log_Map(m,xi(:,j),method);   
        end
    end
    
    if size(xi,3)>1
        s = sum(ss,3);
    else
        s = sum(ss,2);
    end

    m = exp_Map(m,step/N*s,method);   

    err = norm(pm - m);
%     fprintf('mean computation err (iteration %d): %f\n',i,err);
    if err < tol
        break;
    end
    i = i + 1;
    
    if i == 30 % primitive, change to real line search
        step = 0.5;
    end
end

 assert(i <= maxIter);
