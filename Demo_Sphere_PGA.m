%% This demo first shows generating sphere manifold data, and then performs
%% principal geodesic analysis (PGA) to find the mean, eigenvalues and  Eigenvectors
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%  If you use our code, please cite our paper:
%  Zhang, Youshan. "K-means principal geodesic analysis on riemannian manifolds." 
%  In Proceedings of the Future Technologies Conference, pp. 578-589. Springer, Cham, 2019.

%%  Create sphere manifold data
clear;
Dim=3;
numPts=300;% number of points
randn('seed',999);
[~, N]=randOrthoMat(zeros(Dim,Dim));
True_base=N(:,1);
True_w=N(:,2:end);
% True_lambda=diag(rand(size(True_w,2),1));
True_lambda=diag([0.2 0.1]);
True_wLambda = True_w * True_lambda;
sigma=0.1;
x=randn(Dim-1,2*numPts);
acc=1;counter=0;
while(counter<numPts)

 if (norm(True_wLambda * x(:,acc), "fro") < pi/2.0)   
    rgdHat = exp_Map(True_base, True_wLambda * x(:,acc));   
    data(:,acc)=  RandomGenerator(rgdHat, sigma);
    acc=acc+1;
 end
 counter=counter+1;
end
disp('True parameters are:')
True_base
True_lambda
True_w

%% plot the data
r=1;
theta=linspace(0,pi);
phi=linspace(0,2*pi);
[tt,pp]=meshgrid(theta,phi);
[x,y,z] = sph2cart(pp,pi/2-tt,r);
h=figure;
surf(x,y,z)
map=[0.75 0.75 0.75];
axis square
axis off
colormap(gray)
grid off;
shading interp
hold on;
plot3(data(1,:),data(2,:),data(3,:),'.');set(gcf,'color','w');
m=-2.9;
for i=1:1:60

  matToCol = True_w(:,1) * True_lambda(1,1)*m;
  
  tangent= listToTangent(matToCol);

  gtrueth (:,i) = exp_Map(True_base, tangent);
  m=m+0.1;

end
hold on;
h(1)=plot3(gtrueth(1,:),gtrueth(2,:),gtrueth(3,:),'linewidth',5);

pause(3);
%% perform PGA 
[mu, lambda, vector]=PGA(data);


m=-2.9;
for i=1:1:60

  matToCol = vector(:,1) * sqrt(lambda(1,1))*m;
  
  tangent= listToTangent(matToCol);

  estimated (:,i) = exp_Map(mu, tangent);
  m=m+0.1;

end
hold on;
h(2)=plot3(estimated(1,:),estimated(2,:),estimated(3,:),'linewidth',5);

legend(h(1:2),'True geodesic line','Estimated geodesic line')


disp('Estimated parameters are:')
values = sqrt(lambda);
mu
values
vector