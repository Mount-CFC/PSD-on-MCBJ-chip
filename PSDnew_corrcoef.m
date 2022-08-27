clc
clear
close all
tic


fs=10000; % @@@ Sampling Frequency @@@10kHz for MCBJ-Raman
ft=1/fs;
traceL = 500;  %%@@@@@@@@@@ traces length (how many points per trace)@@@@@@@@@@@@@
lowG=-2.8;
highG=-2.5;
lowPSD=-5;
highPSD=-3;
n_bins=50;         %bins of PSD and Gavg
PlotSelect = 0;     %If it is True, polt the raw .TDMS files

[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else
    filename1{1}=filename;
end
num_files = length(filename1)

%%% plot long traces
if PlotSelect == 1
    figure(20)
    if num_files == 1
        test=TDMS_readTDMSFile(filename1{1});
        data_s=test.data{1,4};
        plot(data_s)
        title(filename1{1},'FontSize',5)
        ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
        xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
    else  
        for n = 1:num_files
            test=TDMS_readTDMSFile(filename1{n});
            data_s=test.data{1,4}; 
            subplot(num_files,1,n)
            plot(data_s)
            title(filename1{n},'FontSize',5)
            ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
            xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
        %     clear test data_s

        end
    end
    saveas(gcf,'1_TraceBeforeCut.fig')
end

%% get the entire data
data_s = [];
if num_files == 1
    test=TDMS_readTDMSFile(filename1{1});
    data_s=test.data{1,4};
else  
    for n = 1:num_files
        test=TDMS_readTDMSFile(filename1{n});
        data_temp=test.data{1,4}; %Conductance is the 4th column
        data_s = [data_s,data_temp];
        disp(['loading ' filename1{n} '...']); % Present the file name
    %     clear test data_s

    end
end
%% cut the entire data to single curves

for i=1:(length(data_s)/traceL)
    raw_for_flk{i} = data_s(:, (traceL*(i-1)+1):traceL*i);
end
%% preview 10 of all single curves
figure(1)
for i=1:10
    
    subplot(5,2,i);
    n=unidrnd(length(raw_for_flk));
    plot(raw_for_flk{n});
    title(n)
    ylim([lowG highG])
%     ylim([-3.6 -2.4])
end
saveas(gcf,'2_RawDataForPSD.fig')

%%  FFT and calculate PSD
disp('Data processing...')
for k=1:length(raw_for_flk)
    data_L=2^nextpow2(length(raw_for_flk{k})); %   Determine Number of points (next power of 2)
    yt=transpose(raw_for_flk{k});  % Transpose the column vector to row vector
    yt=power(10,yt.*1);
    Fyt=fft(yt,data_L)/data_L*2; % FFT amplitude
    f_range=fs/data_L*(0:1:data_L-1); % Frequency range
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
    PSD0(k)=trapz(f_range(sp:ep),SQ_Fyt(sp:ep));            % @@@ Trapezoidal integration
    % @@@ correction integration  100 and 1000 are the intergation range
    if f_range(sp)~=100
        PSD0(k)=PSD0(k)+0.5*(((SQ_Fyt(ep+1)-SQ_Fyt(ep))/F*(1000-f_range(ep))+SQ_Fyt(ep))+SQ_Fyt(ep))*(1000-f_range(ep))-0.5*(((SQ_Fyt(sp+1)-SQ_Fyt(sp))/F*(100-f_range(sp))+SQ_Fyt(sp))+SQ_Fyt(sp))*(100-f_range(sp));
    end
    MeanG(k)=log10(mean(yt));    %Gavg
    PSD(k)=log10(PSD0(k)/(10^MeanG(k)));  %PSD normalized by Gavg
  
end


%% plot PSD


figure(2)
Xedges = linspace(lowG,highG,n_bins+1);
Yedges = linspace(lowPSD,highPSD,n_bins+1);   %creat same length with PSD and Gavg

% figure(101)
histcount0 = histcounts2(MeanG,PSD,Xedges,Yedges);
Xaxis = linspace(lowG,highG,n_bins);
Yaxis = linspace(lowPSD,highPSD,n_bins);
histcount = histcount0';

imagesc(Xaxis,Yaxis, histcount);
% imagesc(Xedges,Yedges, histcount);
set(gca,'YDir','normal')
set(gca,'tickdir','out')
load MyColormapRandB.mat;
colormap(mycmap);
colorbar;
xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial') 
ylabel('Noise Power / \itG \rm(log \itG\rm_0)', 'Interpreter','tex','FontSize',20,'FontName','Arial')
set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')

%% two dimensional Gauss fitting
hold on;
figure(2)
[fitresult0, gof] = createFit(Xaxis, Yaxis, histcount);

xRangeStart=lowG;
xRangeEnd=highG;
yRangeStart=lowPSD;
yRangeEnd=highPSD;
    
    
a = fitresult0.a;
b = fitresult0.b;
c = fitresult0.c;
d = fitresult0.d;
e = fitresult0.e;
f = fitresult0.f;

n=1;
meanlogG=mean(MeanG);
meanshift=(1-n)*meanlogG;  %meanshift == 0
[x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
contour(x,y,z,4,'-k','linewidth',1.5);
saveas(gcf,'3_PSD.fig')
    
    


%% determine Gavg^{n}

% lowG=-2.8;
% highG=-2.4;
% lowPSDForN=-5.5;
% highPSDForN=-1;
% Xedges = linspace(lowG,highG,n_bins+1);
% Yedges = linspace(lowPSDForN,highPSDForN,n_bins+1);   
% Xaxis = linspace(lowG,highG,n_bins);
% Yaxis = linspace(lowPSDForN,highPSDForN,n_bins);
PSD_index = 11;
for n = 1.1 : 0.1 : 2
    meanshift=(1-n)*meanlogG;
    for k=1:length(PSD)
        PSD1{PSD_index}(k)=log10(PSD0(k)/(10^MeanG(k))^n);
    end

%     PSD_data1= [transpose(MeanG),transpose(PSD1)];
%     [PSD_2D11,X11,Y11]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
%     [histcount,Xedges,Yedges] = histcounts2(MeanG,PSD1,[n_bins,n_bins]);
    

    Yedges = linspace(lowPSD + meanshift, highPSD + meanshift, n_bins+1);
    Yaxis = linspace(lowPSD + meanshift, highPSD + meanshift, n_bins);
    histcount0 = histcounts2(MeanG,PSD1{PSD_index},Xedges,Yedges);
    histcount = histcount0';
    [fitresult0, gof] = createFit(Xaxis, Yaxis, histcount);

    figure(round(2+10*(n-1)));

    imagesc(Xaxis, Yaxis, histcount);
    set(gca,'YDir','normal')
    set(gca,'tickdir','out')
    load MyColormapRandB.mat;
    colormap(mycmap);
    colorbar;
    xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial') 
    ylabel(['Noise Power / \itG \rm(log \itG\rm_0^{n})  n=' num2str(n)], 'Interpreter','tex','FontSize',20,'FontName','Arial')
    set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')


hold on;

xRangeStart=lowG;
xRangeEnd=highG;
yRangeStart=lowPSD;
yRangeEnd=highPSD;


a = fitresult0.a;
b = fitresult0.b;
c = fitresult0.c;
d = fitresult0.d;
e = fitresult0.e;
f = fitresult0.f;

% meanshift=(1-n)*meanlogG;
[x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
contour(x,y,z,5,'-k','linewidth',1.5);

PSD_index = PSD_index + 1;
disp(num2str(PSD_index));

end


%%  计算PSD和meanG的相关系数，找到最小n
% This is the result determining step!


coe=[];
k=1;
for i = 0.5:0.1:2.5
   
   PSD_temp=log10(PSD0./(10.^MeanG).^i);
   coe_tmp=corrcoef(PSD_temp,MeanG);
   coe(k,1:2)=[i, coe_tmp(1,2)];
   

   k=k+1;
   
end
figure(13)
plot(abs(coe(:, 1)), abs(coe(:, 2)),'-k','Marker','o','MarkerSize',8,'LineStyle','-','LineWidth', 2)
title(['traceL:' num2str(traceL)])
xlabel('\itG\rm_{AVG}^{n}')  % use {} to incorporate all char string
ylabel('Correlation coefficients')
min_corrcoef=coe(abs(coe(:,2))==min(abs(coe(:,2))));

disp(['The scaling exponent zero correlation between the normalized ficker noise power and mean conductance: ',num2str(min_corrcoef)]); 





toc