clc
tic
clear all
close all
load traces_open.mat;


G_Max=-0.3;
G_Min=-5;
blankGmax=-0.3
blankGmin=-1;
n_bins=500;

logGedges=G_Min:(G_Max-G_Min)/n_bins:G_Max
for i=1:length(logG_open)
    h(i,:)=histc(logG_open{i},logGedges);
end
cov_logG=corrcoef(h);
%maxCov=max(max(cov_logG));
%cov_logG=cov_logG/maxCov;
figure(1);
imagesc(logGedges,logGedges, cov_logG);
set(gca,'YDir','normal')
load covmap.mat;
colormap(mycmp);
colorbar;
xlabel('logG (G0)');ylabel('logG (G0)');

toc