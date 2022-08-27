tic
clc
clear all
close all
load raw_for_flk.mat;



fs=50000; % @@@ Sampling Frequency @@@
ft=1/fs;
traceL=10000; %@@@@@@@@@@ traces length (how mant points per trace)@@@@@@@@@@@@@
lowG=-5;
highG=-1;
lowPSD=-5;
highPSD=-2;
n_bins=40;
select=0; %@@@@ 定为0.则把台阶掐头去尾@@@
selectFit=1; %@@@@ 定为0，则开始拟合
cutL=500;


%  Delte bad traces



Lraw = length(raw_for_flk);
z=1;
for i=1:Lraw
    if length(raw_for_flk{z}) <= traceL  % delte traces with lengh shorter than 500 points
        raw_for_flk(z)=[];
        continue;
    end
    z=z+1;
    if i==length(raw_for_flk)
        break;
    end
end

Lraw=length(raw_for_flk);
z=1;
for i=1:Lraw
    if length(raw_for_flk{z}) <= traceL  % delte traces with lengh shorter than 500 points
        raw_for_flk(z)=[];
        continue;
    end
    z=z+1;
    if i==length(raw_for_flk)
        break;
    end
end

t=0;
if t== select
    Lraw=length(raw_for_flk);
    for i=1:Lraw
        raw_for_flk{i}=raw_for_flk{i}(cutL:length(raw_for_flk{i})-cutL);
    end
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FFT and calculate PSD
for k=1:length(raw_for_flk)
    data_L=2^nextpow2(length(raw_for_flk{k})); %   Determine Number of points (next power of 2)
    yt=transpose(raw_for_flk{k});  % Transpose the column vector to row vector
    yt=power(10,yt.*1);
    Fyt=fft(yt,data_L)/data_L*2; % FFT amplitude
    f=fs/data_L*(0:1:data_L-1); % Frequency range
    F=fs/data_L;
    A_Fyt=abs(Fyt);
    SQ_Fyt=A_Fyt.^2; % square
    
    for z=0:data_L              %  intergration strating range
        if fs/data_L*z > 100
            sp=z;
            break
        end
    end
    
    for z=sp:data_L              %  intergration end range
        if fs/data_L*z > 1000
            ep=z;
            break
        end
    end
    %int_range=(100:fs/data_L:1000);
    PSD0(k)=trapz(f(sp:ep),SQ_Fyt(sp:ep));            % @@@ Trapezoidal integration
    % @@@ correction integration  100 and 1000 are the intergation range
    if f(sp)~=100
        PSD0(k)=PSD0(k)+0.5*(((SQ_Fyt(ep+1)-SQ_Fyt(ep))/F*(1000-f(ep))+SQ_Fyt(ep))+SQ_Fyt(ep))*(1000-f(ep))-0.5*(((SQ_Fyt(sp+1)-SQ_Fyt(sp))/F*(100-f(sp))+SQ_Fyt(sp))+SQ_Fyt(sp))*(100-f(sp));
    end
    MeanG(k)=log10(mean(yt));
    PSD(k)=log10(PSD0(k)/(10^MeanG(k)));
    % PSD(k)=PSD(k)/(10^MeanG(k));
    %figure(k);
    % plot(f(1:data_L/2),SQ_Fyt(1:data_L/2));
    % xlim([90,1100]);
end

%maxpsd= max(PSD);
%minpsd= min(PSD);
n=1;
meanlogG=mean(MeanG);
PSD_data= [transpose(MeanG),transpose(PSD)];
[PSD_2D10,X10,Y10]=plotPSD(PSD_data,highPSD,lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
[fitresult0, gof] = createFit(X10, Y10, PSD_2D10);
figure(1);

imagesc(X10,Y10, PSD_2D10);
set(gca,'YDir','normal')
set(gca,'tickdir','out')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (log G0)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (log G0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
    

%% 拟合确定指数值
n=1.1;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD11(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD11)];
[PSD_2D11,X11,Y11]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)

[fitresult0, gof] = createFit(X11, Y11, PSD_2D11);

figure(2);

imagesc(X11,Y11, PSD_2D11);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.1)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@@@@

n=1.2;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD12(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD12)];
[PSD_2D12,X12,Y12]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
[fitresult0, gof] = createFit(X12, Y12, PSD_2D12);
figure(3);

imagesc(X12,Y12, PSD_2D12);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.2)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@
n=1.3;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD13(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD13)];
[PSD_2D13,X13,Y13]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
[fitresult0, gof] = createFit(X13, Y13, PSD_2D13);
figure(4);

imagesc(X13,Y13, PSD_2D13);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.3)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@
n=1.4;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD14(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD14)];
[PSD_2D14,X14,Y14]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(5);
[fitresult0, gof] = createFit(X14, Y14, PSD_2D14);
imagesc(X14,Y14, PSD_2D14);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.4)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@
n=1.5;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD15(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD15)];
[PSD_2D15,X15,Y15]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(6);
[fitresult0, gof] = createFit(X15, Y15, PSD_2D15);
imagesc(X15,Y15, PSD_2D15);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.5)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
n=1.6;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD16(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD16)];
[PSD_2D16,X16,Y16]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(7);
[fitresult0, gof] = createFit(X16, Y16, PSD_2D16);
imagesc(X16,Y16, PSD_2D16);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.6)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@@

n=1.7;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD17(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD17)];
[PSD_2D17,X17,Y17]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(8);
[fitresult0, gof] = createFit(X17, Y17, PSD_2D17);
imagesc(X17,Y17, PSD_2D17);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.7)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@@@@@@
n=1.8;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD18(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD18)];
[PSD_2D18,X18,Y18]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(9);
[fitresult0, gof] = createFit(X18, Y18, PSD_2D18);
imagesc(X18,Y18, PSD_2D18);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.8)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@@@@@@@@
n=1.9;
meanshift=(1-n)*meanlogG
for k=1:length(PSD)
    PSD19(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD19)];
[PSD_2D19,X19,Y19]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(10);
[fitresult0, gof] = createFit(X19, Y19, PSD_2D19);
imagesc(X19,Y19, PSD_2D19);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^1.9)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
%% @@@@@@@@@@@@@@@@@@

n=2;
meanshift=(1-n)*meanlogG;
for k=1:length(PSD)
    PSD2(k)=log10(PSD0(k)/(10^MeanG(k))^n);
end

PSD_data1= [transpose(MeanG),transpose(PSD2)];
[PSD_2D20,X20,Y20]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
figure(11);
[fitresult0, gof] = createFit(X20, Y20, PSD_2D20);
imagesc(X20,Y20, PSD_2D20);
set(gca,'YDir','normal')
load PSD_color.mat;
colormap(PSD_color);
colorbar;
ylabel({'Noise Power/G (logG0^2)'},'Interpreter','latex','FontSize',15)
xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)

hold on;
    xRangeStart=lowG;
    xRangeEnd=highG;
    yRangeStart=lowPSD;
    yRangeEnd=highPSD;
    
    
       a =      fitresult0.a ;
       b =      fitresult0.b;
       c =      fitresult0.c  ;
       d =      fitresult0.d  ;
       e =      fitresult0.e  ;
       f =      fitresult0.f  ;

    meanshift=(1-n)*meanlogG;
    [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
    z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
    contour(x,y,z,5,'-k','linewidth',1.3);
    
%@@@@@@@@@@@@@@@@@ 拟合曲线绘制
% if selectFit==0
%     
%     figure(1);
%     hold on;
%     xRangeStart=lowG;
%     xRangeEnd=highG;
%     yRangeStart=lowPSD;
%     yRangeEnd=highPSD;
%   a =      -6.642 ;
%        b =       5.506  ;
%        c =       -8.99  ;
%        d =      -5.875  ;
%        e =      -2.908  ;
%        f =      -14.65  ;
% 
% 
% 
%     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
%     z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
%     contour(x,y,z,5,'-k','linewidth',1.3);
%     
%     figure(10);
%     hold on;
%     xRangeStart=lowG;
%     xRangeEnd=highG;
%     yRangeStart=lowPSD;
%     yRangeEnd=highPSD;
%        a =      -4.123 ;
%        b =      0.1894 ;
%        c =      -14.67 ;
%        d =        -6.1  ;
%        e =      -2.923  ;
%        f =      -15.14  ;
%  
%     n=2;
%     meanshift=(1-n)*meanlogG;
%     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
%     z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
%     contour(x,y,z,5,'-k','linewidth',1.3);
%     
%     
%     
% end


%@@@@@@@@@@@@@@@@@@@  计算PSD和meanG的相关系数，找到最小的尺度化参数


coe=[];
k=1;
for i = 0.5:0.1:2.5
   
   PSD_temp=log10(PSD0./(10.^MeanG).^i);
   coe_tmp=corrcoef(PSD_temp,MeanG);
   coe(k,1:2)=[i, coe_tmp(1,2)];
   

   k=k+1;
   
end
min_corrcoef=coe(abs(coe(:,2))==min(abs(coe(:,2))))

disp('the scaling exponent zero correlation between the normalized ficker noise power and mean conductance.'); min_corrcoef
toc



    %% plot - fitness @hy
% figure(22)
% subplot(221)
% PSD_G0=log10(PSD0./(10.^MeanG));
% [PSD_2D221,Xedges,Yedges] = histcounts2(MeanG,PSD_G0,[n_bins n_bins],'BinMethod','fd');
% X221=Xedges(1:end-1);
% Y221=Yedges(1:end-1);
% imagesc(X221,Y221, PSD_2D221);
% set(gca,'YDir','normal')
% set(gca,'tickdir','out')
% load PSD_color.mat;
% colormap(PSD_color);
% colorbar;
% ylabel({'Noise Power/G (logG0)'},'Interpreter','latex','FontSize',15)
% xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
% [fitresult1, gof] = createFit(X221, Y221, PSD_2D221);
%  hold on;
%     xRangeStart=lowG;
%     xRangeEnd=highG;
%     yRangeStart=lowPSD;
%     yRangeEnd=highPSD;
%     
%    
%     
%        a =      fitresult1.a ;
%        b =      fitresult1.b;
%        c =      fitresult1.c  ;
%        d =      fitresult1.d  ;
%        e =      fitresult1.e  ;
%        f =      fitresult1.f  ;
%  
% 
% 
% 
%     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,yRangeStart:0.01:yRangeEnd);
%     z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
%     contour(x,y,z,5,'-k','linewidth',1.3);
% 
% 
% n = min_corrcoef;
% 
% meanshift=(1-n)*meanlogG
% % for k=1:length(PSD)
% %     PSD19(k)=log10(PSD0(k)/(10^MeanG(k))^n);
% % end
% PSD_G0=log10(PSD0./(10.^MeanG));
% % PSD_data1= [transpose(MeanG),transpose(PSD19)];
% % [PSD_2D222,X222,Y222]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins);
% % n=1.5
% % meanshift=(1-n)*meanlogG;
% % for k=1:length(PSD)
% %     PSD2(k)=log10(PSD0(k)/(10^MeanG(k))^n);
% % end
% PSD_minum=log10(PSD0./(10.^MeanG).^n);
% 
% % PSD_data1= [transpose(MeanG),transpose(PSD2)];
% 
% % [PSD_min,Xmin,Ymin]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
% 
% [PSD_2D222,Xedges,Yedges] = histcounts2(MeanG,PSD_minum,[n_bins n_bins],'BinMethod','fd');
% X222=Xedges(1:end-1);
% Y222=Yedges(1:end-1);
% figure(22)
% subplot(222)
% % imagesc(Xmin,Ymin, PSD_min);
% imagesc(X222,Y222, PSD_2D222);
% set(gca,'YDir','normal')
% set(gca,'tickdir','out')
% load PSD_color.mat;
% colormap(PSD_color);
% colorbar;
% str=['Noise Power/G (logG0^' num2str(n) ')'];
% ylabel({str},'Interpreter','latex','FontSize',15)
% xlabel({'Conductance (logG0)'},'Interpreter','latex','FontSize',15)
% [fitresult2, gof] = createFit(X222, Y222, PSD_2D222);
% hold on;
%     xRangeStart=lowG;
%     xRangeEnd=highG;
%     yRangeStart=lowPSD;
%     yRangeEnd=highPSD;
%     
%     
%        a =      fitresult2.a ;
%        b =      fitresult2.b;
%        c =      fitresult2.c  ;
%        d =      fitresult2.d  ;
%        e =      fitresult2.e  ;
%        f =      fitresult2.f  ;
% 
%     meanshift=(1-n)*meanlogG;
%     [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
%     z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
%     contour(x,y,z,5,'-k','linewidth',1.3);
% toc