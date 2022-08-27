% select one sigle data to obtain the hist plot
clc
clear 
tic


[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else filename1{1}=filename;
end


test=TDMS_readTDMSFile(filename1{1});
data_s=test.data{1,4}; %Conductance is in the 4th column
% histogram(data_s, 1000)  %all histogram
lower_time = 589000
upper_time = 633000
figure
subplot(121)
plot(data_s(:,lower_time:upper_time))
title(filename1{1})
ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',12)
% title({['abc','L Range:',num2str(a),'(nm)','~~',num2str(b),'(cm)'];['B Range:',num2str(c),'(cm)','~~',num2str(d),'(cm)']})
% https://www.cnblogs.com/AI-Algorithms/p/3731232.html  
% Dont forget: num2str
xlabel({['Sampling points'];['From ' num2str(lower_time)  ' to '  num2str(upper_time)]},'Interpreter','tex','FontSize',12)

subplot(122)
histogram(data_s(:,lower_time:upper_time), 1000)
xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',12)
ylabel({'Counts'},'Interpreter','tex','FontSize',12)
toc