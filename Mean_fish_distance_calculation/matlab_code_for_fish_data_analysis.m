close all
clear                                                                      %%%% This clears all the variables in the workspace
clc                                                                        %%%% This clears all the texts from the command window

addpath('/Users/pradeepkumar/Box Sync/PhD work/matlab/manuscript/zscore');     %%%% path where this file is present
cd('/Users/pradeepkumar/Box Sync/PhD work/matlab/manuscript/zscore');          %%%% path where the data file is present

%%%% loading the data
[num,txt,data] = xlsread('fish.xlsx');
txt_data=txt(2:end,1);

%%%% finding the unique genomic coordinates in the dataset
unique_genomic_coordinates=unique(txt_data);

data_without_header=data(2:end,:);
size_all_data=size(data_without_header);
size_of_unique_data=size(unique_genomic_coordinates);

%%%% saving the whole data into cells of repeating genomic data
for i=1:(size_all_data(1)/size_of_unique_data(1))
    data_info_cell{1,i}=data_without_header(1+(size_of_unique_data(1)*(i-1)):size_of_unique_data(1)+(size_of_unique_data(1)*(i-1)),:);
end

data=data;

%%%% separating all the datasets now into separate matrices
for mn=1:size(data_info_cell,2)
    genomic_coordinate_data(:,mn)=data_info_cell{1,mn}(:,1);
    distance_to_lamina(:,mn)=data_info_cell{1,mn}(:,2);
    distance_to_nucleoli(:,mn)=data_info_cell{1,mn}(:,3);
    distance_to_speckles(:,mn)=data_info_cell{1,mn}(:,4);
    test_genomic_data_arrays(:,mn)=strcmp(genomic_coordinate_data(:,1),genomic_coordinate_data(:,mn));
end

data=data;
%%
%%%% finiding mean distance from lamina. I needed to convert nan to
%%%% character and then calculate nanmean. 
distance_to_lamina(cellfun(@ischar,distance_to_lamina)) = {nan};
distance_to_lamina  = cell2mat(distance_to_lamina);
mean_distance_to_lamina=nanmean(distance_to_lamina,2);
std_distance_to_lamina=nanstd(distance_to_lamina,0,2);

for pq=1:size(distance_to_lamina,1)
    zscore_lamina(pq,:)=(distance_to_lamina(pq,:)-mean_distance_to_lamina(pq,1))./std_distance_to_lamina(pq,1);
end

mean_zscore_distance_to_lamina=nanmean(zscore_lamina,2);

%%%% finiding mean distance from nucleoli. I needed to convert nan to
%%%% character and then calculate nanmean. 
distance_to_nucleoli(cellfun(@ischar,distance_to_nucleoli)) = {nan};
distance_to_nucleoli  = cell2mat(distance_to_nucleoli);
mean_distance_to_nucleoli=nanmean(distance_to_nucleoli,2);

std_distance_to_nucleoli=nanstd(distance_to_nucleoli,0,2);

for pq=1:size(distance_to_nucleoli,1)
    zscore_nucleoli(pq,:)=(distance_to_nucleoli(pq,:)-mean_distance_to_nucleoli(pq,1))./std_distance_to_nucleoli(pq,1);
end

mean_zscore_distance_to_nucleoli=nanmean(zscore_nucleoli,2);


%%%% finiding mean distance from speckles. I needed to convert nan to
%%%% character and then calculate nanmean. 
distance_to_speckles(cellfun(@ischar,distance_to_speckles)) = {nan};
distance_to_speckles  = cell2mat(distance_to_speckles);
mean_distance_to_speckles=nanmean(distance_to_speckles,2);

%%%% making a final matrix
final_data=[genomic_coordinate_data(:,1) num2cell(mean_distance_to_lamina)...
    num2cell(mean_distance_to_nucleoli) num2cell(mean_distance_to_speckles)];

%%%% adding header row to the matrix
final_data_with_header=[data(1,:)
    final_data];

%%%%% writing the matrix final_data into an excel file
writecell(final_data_with_header,'fish_data_with_mean_upload.xlsx')