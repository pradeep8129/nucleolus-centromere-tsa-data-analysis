%%% Scatter plot code for Pradeep
clear
clc
close all

addpath(%'path where this MATLAB file is present');    
cd(%'path where the data file to be analyzed is present');    

[num,txt,raw] = xlsread('fish_data_with_mean.xlsx'); 
%%
original_text=txt(1,:);
txt=txt(2:end,:);
index_of_interest=strcmp(txt(:,1),'chr4'); % put chr name here 
data_of_interest=num([index_of_interest],:);

scatter_marker_size=50;
x_column=3;
y_column=4;
coordinate_of_interest=16; % The coordinate of interest is the last fish probe for p arm. For example for chromosome 4, first 16 probes are for p arm and after that the probes are q-arm 

figure1 = figure;
figure1.Position = [500 400 800 600];  % 560   420
scatter(num(:,x_column),num(:,y_column),scatter_marker_size,'MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerFaceAlpha',0.2, 'MarkerEdgeAlpha',0.2)
hold on;
scatter(data_of_interest(1:coordinate_of_interest,x_column),data_of_interest(1:coordinate_of_interest,y_column),scatter_marker_size,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold on;
scatter(data_of_interest(coordinate_of_interest+1:end,x_column),data_of_interest(coordinate_of_interest+1:end,y_column),scatter_marker_size,'MarkerEdgeColor','b','MarkerFaceColor','b')
hold on;

xlabel(original_text(x_column+1),'fontsize',20,'Interpreter','none');
ylabel(original_text(y_column+1),'fontsize',20,'Interpreter','none');               %%% title
box on;                                                 
set(gca,'fontsize',30)                                  

ylim([1.6 2.7])
xlim([0.32 0.75])

legend('','chr 4 p-arm','chr 4 q-arm')

