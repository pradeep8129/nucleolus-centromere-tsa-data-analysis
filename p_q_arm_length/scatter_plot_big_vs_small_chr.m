clear
clc
close all

addpath(%'path where this MATLAB file is present');     
cd(%'path where the data file to be analyzed is present');     

[num,txt,raw] = xlsread('fish_data_with_mean.xlsx'); 
%%
original_text=txt(1,:);
txt=txt(2:end,:);
index_of_interest_1=strcmp(txt(:,1),'chr1');
data_of_interest_1=num([index_of_interest_1],:);

index_of_interest_2=strcmp(txt(:,1),'chr2');
data_of_interest_2=num([index_of_interest_2],:);

index_of_interest_3=strcmp(txt(:,1),'chr3');
data_of_interest_3=num([index_of_interest_3],:);

index_of_interest_4=strcmp(txt(:,1),'chr18');
data_of_interest_4=num([index_of_interest_4],:);

index_of_interest_5=strcmp(txt(:,1),'chr19');
data_of_interest_5=num([index_of_interest_5],:);

index_of_interest_6=strcmp(txt(:,1),'chr20');
data_of_interest_6=num([index_of_interest_6],:);

scatter_marker_size=50;
x_column=3;
y_column=4;

figure1 = figure;
figure1.Position = [500 400 800 600];  % 560   420
scatter(num(:,x_column),num(:,y_column),scatter_marker_size,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerFaceAlpha',0.3, 'MarkerEdgeAlpha',0.3)
hold on;
scatter(data_of_interest_1(:,x_column),data_of_interest_1(:,y_column),scatter_marker_size,'MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD')
hold on;
scatter(data_of_interest_2(:,x_column),data_of_interest_2(:,y_column),scatter_marker_size,'MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
hold on;
scatter(data_of_interest_3(:,x_column),data_of_interest_3(:,y_column),scatter_marker_size,'MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120')
hold on;
scatter(data_of_interest_4(:,x_column),data_of_interest_4(:,y_column),scatter_marker_size,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30')
hold on;
scatter(data_of_interest_5(:,x_column),data_of_interest_5(:,y_column),scatter_marker_size,'MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F')
hold on;
scatter(data_of_interest_6(:,x_column),data_of_interest_6(:,y_column),scatter_marker_size,'MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE')


xlabel(original_text(x_column+1),'fontsize',20,'Interpreter','none');
ylabel(original_text(y_column+1),'fontsize',20,'Interpreter','none');               %%% title
box on;                                                 
set(gca,'fontsize',30)                                  
ylim([0.5 2.7])
xlim([0.3 0.8])

legend('','chr1','chr2','chr3','chr18','chr19','chr20')

