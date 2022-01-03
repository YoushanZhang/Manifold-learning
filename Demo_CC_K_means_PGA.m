%% K means PGA on two clusters corpus callosum data on Kendall's shape manifold
%% Created by Youshan Zhang
%  Last modified: 09/30/2021
%  If you have any questions, please contact me at zys618@hotmail.com.
%%  If you use our code, please cite our paper:
%  Zhang, Youshan. "K-means principal geodesic analysis on riemannian manifolds." 
%  In Proceedings of the Future Technologies Conference, pp. 578-589. Springer, Cham, 2019.

%%
clear;
clc;
close all;
load('ccdata.mat');
mix=2;

[pi, mu, vector, lambad, sigma2,clusters] = k_means_PGA(data, mix,'Kendall', 1);
 m=-1:0.1:1;
cmap = cool(size(m,2));   
for k=1:mix
  figure,  
for i=1:size(m,2)
 y=exp_Map(mu(:,:,k),reshape(vector(:,1,k),[2,size(data,2)])*...
           sqrt(lambad(1,1,k))*m(i), 'Kendall');
 plot(y(1,:),y(2,:),'Color', cmap(i,:))
 hold on;
end
set(gcf,'color','w')
set(gca,'ydir','reverse');
axis off; 
end

colormap(cmap);
set(gcf,'color','w');
colorbar('Direction','reverse','Ticks',[0,0.5,1],...
         'TickLabels',{'-\lambda','0','\lambda'})