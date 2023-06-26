%% group ROI wrapper
% this script contains the group analyses for all the ROI results 
% 1. Manmade/Natural decoding 
% 2. Distance-to-hyperplane correlation with categorization RTs
% 3. Distance-to-hyperplane correlation with distraction RTs

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

% setup spm
try
    spm;
    close all
catch
    spm_path = '/scratch/singej96/dfg_projekt/WP1/analysis_tools/spm12';
    addpath(spm_path);
end


% set plot defaults 

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica') 

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';
% get fmri subnames 

fmri_subs = dir(fullfile(results_path,'*sub*')); 

fmri_subs = {fmri_subs.name}';

% specify excluded subjects
excluded_subjects = {'sub12'}; 

%% manmade/natural decoding

decoding_roi = [];

%specify results name
res_name = 'manmade_natural';

for i_sub = 1:length(fmri_subs)
    
    sub_id = fmri_subs{i_sub};
    
    if any(ismember(excluded_subjects, sub_id)), continue, end 
    
    if ~isdir(fullfile(results_path,sub_id,'decoding',res_name,'roi')); fprintf('Not %i\n',i_sub);  continue; end; 
    
    results_dir =  fullfile(results_path,sub_id,'decoding',res_name,'roi');
    load(fullfile(results_dir, 'res_accuracy_minus_chance.mat'));
    decoding_roi = [results.accuracy_minus_chance.output; decoding_roi];
end 

fprintf('Mean decoding accuracies over all subjects EVC: %2f, PPA: %2f, LOC: %2f\n', mean(decoding_roi(:,1)),mean(decoding_roi(:,2)),mean(decoding_roi(:,3))); 

% plot 

roi_names = {'EVC'; 'LOC';'PPA'};

%cmap = colormap('inferno');
cmap = colormap('redblueTecplot');
close all

all_accs = mean(decoding_roi(:,1:3))';
decoding_se = [std(decoding_roi(:,1))/sqrt(length(decoding_roi)),...
            std(decoding_roi(:,2))/sqrt(length(decoding_roi)),...
            std(decoding_roi(:,3))/sqrt(length(decoding_roi))]';

figure
h = bar(all_accs, 'grouped','FaceColor', 'flat');
h.CData= [[0 0 0];cmap(256,:);cmap(200,:)];
xticklabels([roi_names])
yticks([0:5:35])
yticklabels([50:5:85])
xlabel('ROI')
ylabel('Decoding Accuracy (%)')
title('Manmade vs. Natural Decoding')

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchxange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i), decoding_se(:,i), 'k', 'linestyle', 'none', 'linewidth', 2);
end

% loop through each roi
for roi = 1:3
    % plot individual data points for each roi
    scatter(x(roi) * ones(length(decoding_roi), 1), decoding_roi(:,roi), 10, [0.6 0.6 0.6], 'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
end
%legend({'Photos'; 'Drawings'; 'Sketches'} ,'Location','northeast')

print(fullfile(out_dir, ['manmade_natural_decoding_ROI.svg']), ...
             '-dsvg', '-r600')
   
         
%% Distance-to-hyperplane correlation 

dth_corr = [];       

% load behavior 
load(fullfile(behav_path,'RT_all_subjects_5_35_categorization.mat'), 'RTs')

mean_RTs = nanmean(RTs,1);%.*[ones(1,30), ones(1,30)*-1]; 

load(fullfile(behav_path,'RT_all_subjects_5_35_fixation.mat'), 'RTs')

distraction_RTs = nanmean(RTs,1);%.*[ones(1,30), ones(1,30)*-1]; 

%specify results name
res_name = 'manmade_natural';

for i_sub = 1:length(fmri_subs)
    
    sub_id = fmri_subs{i_sub};
    if any(ismember(excluded_subjects, sub_id)), continue, end 
    
    if ~isdir(fullfile(results_path,sub_id,'decoding',res_name,'roi')); fprintf('Not %i\n',i_sub);  continue; end; 
    
    results_dir =  fullfile(results_path,sub_id,'decoding',res_name,'roi');
    load(fullfile(results_dir, 'res_mean_decision_values.mat'));
    
    for i = 1:length(results.mean_decision_values.output)
    
    these_dec_vals = results.mean_decision_values.output{i}; %reshape(results.mean_decision_values.output{i},2,60);
    dec_vals(i_sub,i,:) = these_dec_vals; %store the decision values for later
    if length(these_dec_vals) > 60
        these_dec_vals = mean(reshape(these_dec_vals,length(these_dec_vals)/60,60))';
    end 
    dth_corr(i_sub,i) = corr(these_dec_vals,mean_RTs', 'Type','Pearson');
    dth_corr_distraction(i_sub,i) = corr(these_dec_vals,distraction_RTs', 'Type','Pearson');
    end 


end 

% remove the row of the excluded subjects 
dec_vals(find(dec_vals(:,1) ==0),:,:) = []; 
dth_corr(find(dth_corr(:,1) ==0),:,:) = []; 
dth_corr_distraction(find(dth_corr_distraction(:,1)==0),:,:) = []; 


fprintf('\nMean distance to hyperplane correlation over all subjects EVC: %2f, LOC: %2f, PPA: %2f\n',mean(dth_corr)); 
fprintf('\nMean distance to hyperplane correlation for distraction over all subjects EVC: %2f, LOC: %2f, PPA: %2f\n', mean(dth_corr_distraction)); 

% plot 

roi_names = {'EVC'; 'LOC';'PPA'};%;'OPA';'RSC'};

colorspace = linspace(20,size(cmap,1),length(roi_names)); 
close all

all_accs = mean(dth_corr(:,1:3))';
decoding_se = [];
for roi = 1:length(roi_names)
    decoding_se(roi,:) = [std(dth_corr(:,roi))/sqrt(length(dth_corr))];
end

figure
h = bar(all_accs, 'grouped','FaceColor', 'flat');
h.CData= [[0 0 0];cmap(256,:);cmap(200,:)];
xticklabels([roi_names])
%yticks([0:-0.02:-0.12])
xlabel('ROI')
ylabel('Pearson R')
title('Distance To Hyperplane Correlation')

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i), decoding_se(:,i), 'k', 'linestyle', 'none','linewidth',2);
end


% loop through each roi
for roi = 1:3
    % plot individual data points for each roi
    scatter(x(roi) * ones(length(dth_corr), 1), dth_corr(:,roi), 10, [0.6 0.6 0.6], 'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
end

print(fullfile(out_dir,['dth_corr_ROI.svg']), ...
             '-dsvg', '-r600')
 
% do the same for the distraction DTH-correlations
all_accs = mean(dth_corr_distraction(:,1:3))';
decoding_se = [];
for roi = 1:length(roi_names)
    decoding_se(roi,:) = [std(dth_corr_distraction(:,roi))/sqrt(length(dth_corr_distraction))];
end
           
figure
h = bar(all_accs, 'grouped','FaceColor', 'flat');
h.CData= [[0 0 0];cmap(256,:);cmap(200,:)];
xticklabels([roi_names])
%yticks([0:-0.02:-0.12])
xlabel('ROI')
ylabel('Pearson R')
title('Distance To Hyperplane Correlation Distraction')

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(all_accs, 1);
nbars = size(all_accs, 2);
% Calculate the width for each bar group

groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, all_accs(:,i), decoding_se(:,i), 'k', 'linestyle', 'none','linewidth',2);
end


% loop through each roi
for roi = 1:3
    % plot individual data points for each roi
    scatter(x(roi) * ones(length(dth_corr), 1), dth_corr_distraction(:,roi), 10, [0.6 0.6 0.6], 'filled','MarkerEdgeColor',[0.6 0.6 0.6]);
end

print(fullfile(out_dir,['dth_corr_distraction_ROI.svg']), ...
             '-dsvg', '-r600')
 
%% compute statistics

% set rng 
rng(96)

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

sig_decoding_EVC = permutation_1sample_alld (decoding_roi(:,1), nperm, cluster_th, significance_th, tail);

sig_decoding_LOC = permutation_1sample_alld (decoding_roi(:,2), nperm, cluster_th, significance_th, tail);

sig_decoding_PPA = permutation_1sample_alld (decoding_roi(:,3), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_decoding] = fdr_bh([sig_decoding_EVC sig_decoding_LOC sig_decoding_PPA]);

tail = 'both';

sig_dth_corr_EVC = permutation_1sample_alld(dth_corr(:,1), nperm, cluster_th, significance_th, tail);

sig_dth_corr_LOC = permutation_1sample_alld(dth_corr(:,2), nperm, cluster_th, significance_th, tail);

sig_dth_corr_PPA = permutation_1sample_alld(dth_corr(:,3), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_dth_corr] = fdr_bh([sig_dth_corr_EVC sig_dth_corr_LOC sig_dth_corr_PPA]);

sig_dth_corr_distraction_EVC = permutation_1sample_alld(dth_corr_distraction(:,1), nperm, cluster_th, significance_th, tail);

sig_dth_corr_distraction_LOC = permutation_1sample_alld(dth_corr_distraction(:,2), nperm, cluster_th, significance_th, tail);

sig_dth_corr_distraction_PPA = permutation_1sample_alld(dth_corr_distraction(:,3), nperm, cluster_th, significance_th, tail);

[~,~,~,adj_p_dth_corr_distraction] = fdr_bh([sig_dth_corr_distraction_EVC sig_dth_corr_distraction_LOC sig_dth_corr_distraction_PPA]);

