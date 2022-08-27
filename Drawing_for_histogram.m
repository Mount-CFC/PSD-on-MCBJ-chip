% select one sigle data to obtain the hist plot
clc
clear 
tic


[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else 
    filename1{1}=filename;
end


test=TDMS_readTDMSFile(filename1{1});
data_s=test.data{1,4}; %Conductance is in the 4th column
% histogram(data_s, 1000)  %all histogram

%%
% lower_time = 625000
% upper_time = 626100

lower_time = 1;
upper_time = length(data_s);

conductance = data_s(:,lower_time:upper_time);
time = (1:length(conductance))*1E-4;
figure(1)
plot(time, conductance, 'Color', 'k','LineWidth',1)
% title(filename1{1})
ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial')
% title({['abc','L Range:',num2str(a),'(nm)','~~',num2str(b),'(cm)'];['B Range:',num2str(c),'(cm)','~~',num2str(d),'(cm)']})
% https://www.cnblogs.com/AI-Algorithms/p/3731232.html  
% Dont forget: num2str
% xlabel({['Sampling points'];['From ' num2str(lower_time)  ' to '  num2str(upper_time)]},'Interpreter','tex','FontSize',12)
xlabel(['Time / s'],'Interpreter','tex','FontSize',20,'FontName','Arial')
xlim tight
set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')

% set(gca,'XTick',[])
% box off
% ax2 = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% set(ax2,'YTick', [],'LineWidth',1.5);
% set(ax2,'XTick', [],'LineWidth',1.5);
% box on
% set(gcf,'unit','centimeters','position',[30 5 10 5])

%%
figure(2)
% histogram(conductance, 100,'FaceColor','#D95319','LineStyle','none','Orientation','horizontal');
histogram(conductance, 100,'FaceColor','#D95319','LineStyle','none');

xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial')
ylabel({'Counts'},'Interpreter','tex','FontSize',20,'FontName','Arial')


set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')

toc