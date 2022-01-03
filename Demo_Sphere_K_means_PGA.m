%% K means PGA on three clusters sphere manifold
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%  If you use our code, please cite our paper:
%  Zhang, Youshan. "K-means principal geodesic analysis on riemannian manifolds." 
%  In Proceedings of the Future Technologies Conference, pp. 578-589. Springer, Cham, 2019.
clear
%% Plot the data
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
load data_new.mat
data1=data(:,1:200);data2=data(:,201:489);data3=data(:,490:764);
plot3(data1(1,:),data1(2,:),data1(3,:),'g.'); hold on;plot3(data2(1,:),data2(2,:),data2(3,:),'r.');
hold on; plot3(data3(1,:),data3(2,:),data3(3,:),'k.');


 
%% K means PGA 
load tgeodisc.mat 
hold on;
h(1)=plot3(tgeodisc(1,:),tgeodisc(2,:),tgeodisc(3,:),'y','linewidth',5); % true geodesic line



[pi, mu, vector, lambda, sigma2,clusters] = k_means_PGA(data, 3);
%% Plot estimated line    
for k=1:1:3
m=-1.9;
for i=1:1:40

  matToCol = vector(:,1,k) * sqrt(lambda(1,1,k))*m;
  
  tangent= listToTangent(matToCol);

  estimated (:,i) = exp_Map(mu(:,k), tangent);
  m=m+0.1;

end
hold on;
h(2)=plot3(estimated(1,:),estimated(2,:),estimated(3,:),'b','linewidth',5);
end
set(gcf,'color','w')
legend(h(1:2),'True geodesic line','Estimated geodesic line')