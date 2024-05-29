clear                                               %%% deletes all the past variable in workspace
clc                                                 %%% cleans the text in command window
close all                                           %%% closes all the previous plots

addpath(%'path where this MATLAB file is present');    
cd(%'path where the data file to be analyzed is present'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%% loading the dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%
[num,txt,raw] = xlsread('K562_combined data frame_short_name.xlsx'); 
%%% the header information in the excel file is in variable 'txt'
%%% numerical values in the excel file is in variable 'num'
%%% the whole data in the excel file is in variable 'raw'

%%
all_genetic_locations=raw(2:end,8);
[Grp_Counts,Grp_unique] = groupcounts(all_genetic_locations); 
%%% Grp_counts is the number of time Grp_unique variable appears in an array

%%% sort them in ascending order
[~,sorted_index]=sort(Grp_Counts,'ascend');
sorted_Grp_unique=Grp_unique(sorted_index); %%% I am sorting the genetic location based on the index of Grp_counts

%% using for loop to collect all the value of score for each genetic location
for i=1:size(sorted_Grp_unique,1)
    dummy_indices_of_genetic_location=[];
    dummy_indices_of_genetic_location=strcmp(raw(:,8),sorted_Grp_unique{i,:});
    
    unique_genetic_location_names{:,i}=sorted_Grp_unique{i,:}; %% genetic location names are stored here
    
    dummy_variable=raw(dummy_indices_of_genetic_location,7);
    double_dummy_variable=cell2mat(dummy_variable);
    data_unique_genetic_location{:,i}=double_dummy_variable; %% corresponding scores for each genetic location is stored here
end

%% plotting box plot

x_data=[];
y_data=[];

for k=1:size(sorted_Grp_unique,1)

x_data = [x_data
    k*ones(size(data_unique_genetic_location{1,k},1),1)];

y_data=[y_data
    data_unique_genetic_location{:,k}];

end
%%
f = figure;
f.Position = [100 100 1200 800]; %% size of figure; length=1200, width=800
bh=boxplot(y_data,x_data, 'symbol', ''); %%%  'symbol', '' assign no shape to outliers
% ylabel('MKI67IP TSA-Seq Score')
xticklabels(unique_genetic_location_names)
xtickangle(90)     %% roatating the xtick labels by 90 degree

ylim([-1.5 2.7])

colors =jet(size(unique_genetic_location_names,2));
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',0.8,'Edgecolor','k'); %% if you want darker facecolor then instead of 0.5 use 1
    lines_1 = findobj(gcf, 'type', 'line', 'Tag', 'Box');
    set(lines_1, 'Color', 'k');    %% Median line is given by red color
end
set(gca,'fontsize',20)                                  %%% works like a master command: changes the fontsize for all the objects in the plot
set(bh,'LineWidth', 2);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');%% Median line is given by
