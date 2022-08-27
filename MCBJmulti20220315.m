%生成多张的logG-t图用于粗筛
%一张图生成多张子图的logG-t图用于粗筛

clc
clear 
close all
tic



[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else filename1{1}=filename;
end

num_files = length(filename1)

if num_files == 1
    plot(TDMS_readTDMSFile(filename1{1}).data{1,4})
    title(filename1{1},'FontSize',5)
    ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
    xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
else  
    for n = 1:num_files
        test=TDMS_readTDMSFile(filename1{n});
        data_s=test.data{1,4}; %第一行第4列，提取Conductance
        subplot(num_files,1,n)
        plot(data_s)
        title(filename1{n},'FontSize',5)
        ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
        xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
    %     clear test data_s

    end
end



toc