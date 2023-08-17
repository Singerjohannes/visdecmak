%% group behavioral wrapper
% this script contains the group analyses for all the behavioral results
% computes the correlation between the mean RTs in the categorization and
% distraction task and plots the results 

clear all 
%clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');
% create figure path
if ~isdir(figure_path); mkdir(figure_path); end 
% specify path where first level results are stored
results_path = fullfile(path, 'results'); 
% specify path where neural network modelling results are stored 
modelling_path = fullfile(path,'modelling');
% specify path where behavioral results are stored
behav_path = fullfile(path, 'behav');
% specify where group level results should be stored
out_dir = fullfile(results_path,'group');
if ~isdir(out_dir), mkdir(out_dir), end 

% add stats and utils functions 
addpath(fullfile(pwd,'utils'));
addpath(genpath(fullfile(pwd,'stats')));


% set plot defaults 

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica') 

% load behavior 
load(fullfile(behav_path,'RT_all_subjects_5_35_categorization.mat'), 'RTs')

mean_RTs = nanmean(RTs,1);%.*[ones(1,30), ones(1,30)*-1]; 

load(fullfile(behav_path,'RT_all_subjects_5_35_fixation.mat'), 'RTs')

distraction_RTs = nanmean(RTs,1);%.*[ones(1,30), ones(1,30)*-1]; 

%% compute correlation and plot 

[r,p_val] = corr(mean_RTs',distraction_RTs'); 

figure 

scatter(mean_RTs,distraction_RTs,"filled",'black')
lsline
ylim([0.43 0.47])
yticks([0.43:0.01:0.470])
yticklabels([430:10:470])
ylabel('RT - distraction task (ms)')
xlim([0.47 0.51])
xticks([0.47:0.01:0.51])
xticklabels([470:10:510])
xlabel('RT - categorization task (ms)')
hold on
txt = ['R = ',num2str(r,'%.3f'),'; p = ',num2str(p_val,'%.3f')];
text(0.495,0.457,txt,'FontSize', 14)

print(fullfile(results_path,['categorization_distraction_RT_corr.svg']), ...
             '-dsvg', '-r600')
