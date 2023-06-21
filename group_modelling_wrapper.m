%% group searchlight wrapper
% this script contains the group analyses for the neural network modelling results 

clear all 
clc

%setup paths 

path = fileparts(pwd);
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

%% load fmri and behavioral results 

% get fmri subnames 

fmri_subs = dir(fullfile(results_path,'*sub*'));
fmri_subs = {fmri_subs.name}';

% specify excluded subjects
excluded_subjects = {'sub12'};%{'sub08';'sub14';'sub15';'sub23';'sub29'}; 

dth_corr = [];

% load behavior 
load(fullfile(behav_path,'RT_all_subjects_5_35_categorization.mat'), 'RTs')
mean_RTs = nanmean(RTs,1);

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
    end 
end

dec_vals(find(dec_vals(:,1,1) ==0),:,:) = [];

%% load the model distances 

% load ResNet18 distances 
dist_files = dir(fullfile(modelling_path,'*.mat'));
dist_filenames = {dist_files.name}'; 
filenames_info = regexpi(dist_filenames,'^resnet18_places.*distances.*'); 

resnet18_distances = [];
for idx = 1:length(dist_filenames)
    if filenames_info{idx}
        results = load(fullfile(modelling_path,dist_filenames{idx}));
        resnet18_distances= cat(2,resnet18_distances,abs(results.dec_vals'));
    end 
end 

order_idx = [6 2 3 4 5 1]; % ordered values from low to high level features
resnet18_distances =resnet18_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

% load ResNet50 distances 
dist_files = dir(fullfile(modelling_path,'*.mat'));
dist_filenames = {dist_files.name}'; 
filenames_info = regexpi(dist_filenames,'^resnet50_places.*distances.*'); 

resnet50_distances = [];
for idx = 1:length(dist_filenames)
    if filenames_info{idx}
        results = load(fullfile(modelling_path,dist_filenames{idx}));
        resnet50_distances= cat(2,resnet50_distances,abs(results.dec_vals'));
    end 
end 

order_idx = [6 2 3 4 5 1]; % ordered values from low to high level features
resnet50_distances =resnet50_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

% load AlexNet distances 
dist_files = dir(fullfile(modelling_path,'*.mat'));
dist_filenames = {dist_files.name}'; 
filenames_info = regexpi(dist_filenames,'^alexnet_places.*distances.*'); 

alexnet_distances = [];
for idx = 1:length(dist_filenames)
    
    if filenames_info{idx}
    results = load(fullfile(modelling_path,dist_filenames{idx}));
    alexnet_distances = cat(2,alexnet_distances,abs(results.dec_vals'));
    end 
end 

order_idx = [2 3 5 6 4 1]; % ordered values from low to high level features
alexnet_distances =alexnet_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

% load densenet distances

dist_files = dir(fullfile(modelling_path,'*.mat'));
dist_filenames = {dist_files.name}'; 
filenames_info = regexpi(dist_filenames,'^densenet161_places.*'); 

vgg_distances = [];
for idx = 1:length(dist_filenames)
    
    if filenames_info{idx}
    results = load(fullfile(modelling_path,dist_filenames{idx}));
    vgg_distances = cat(2,vgg_distances,abs(results.dec_vals'));
    end 
end 

order_idx = [6 1 2 3 4 5]; % ordered values from low to high level features
vgg_distances =vgg_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

%% shared variance for each model independently - ResNet18

addpath('/Users/johannessinger/scratch/dfg_projekt/WP1/analysis/utils')

disp('Running commonality analysis...')

shared_var = [];
for sub = 1:size(dec_vals,1)
for j_roi = 1:3
    for model_idx = 1:size(resnet18_distances,2)
        
    %for sub = 1:30 
    xMRI = squeeze(dec_vals(sub,j_roi,:)); 
    xmodel = resnet18_distances(:,model_idx);
    
    y = mean_RTs';
    
    % and now again with another calculation
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    
    shared_var(sub,j_roi,model_idx,1) = r2_mri(1)+r2_model(1)-r2_full(1); 
   % end 
end
end 
end
disp('done.')


%% shared variance for each model independently - ResNet50

addpath('/Users/johannessinger/scratch/dfg_projekt/WP1/analysis/utils')

disp('Running commonality analysis...')


for sub = 1:size(dec_vals,1)
for j_roi = 1:3
    for model_idx = 1:size(resnet50_distances,2)
        
    %for sub = 1:30 
    xMRI = squeeze(dec_vals(sub,j_roi,:)); 
    xmodel = resnet50_distances(:,model_idx);
    
    y = mean_RTs';
    
    % and now again with another calculation
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    
    shared_var(sub,j_roi,model_idx,2) = r2_mri(1)+r2_model(1)-r2_full(1); 
   % end 
end
end 
end
disp('done.')

%% shared variance for each model independently - AlexNet

addpath('/Users/johannessinger/scratch/dfg_projekt/WP1/analysis/utils')

disp('Running commonality analysis...')

for sub = 1:size(dec_vals,1)
for j_roi = 1:3
    for model_idx = 1:size(alexnet_distances,2)
        
    %for sub = 1:30 
    xMRI = squeeze(dec_vals(sub,j_roi,:)); 
    xmodel = alexnet_distances(:,model_idx);
    
    y = mean_RTs';
    
    % and now again with another calculation
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    
    shared_var(sub,j_roi,model_idx,3) = r2_mri(1)+r2_model(1)-r2_full(1); 
   % end 
end
end 
end
disp('done.')

%% shared variance for each model independently - DenseNet

addpath('/Users/johannessinger/scratch/dfg_projekt/WP1/analysis/utils')

disp('Running commonality analysis...')

for sub = 1:size(dec_vals,1)
for j_roi = 1:3
    for model_idx = 1:size(vgg_distances,2)
        
    %for sub = 1:30 
    xMRI = squeeze(dec_vals(sub,j_roi,:)); 
    xmodel = vgg_distances(:,model_idx);
    
    y = mean_RTs';
    
    % and now again with another calculation
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    
    shared_var(sub,j_roi,model_idx,4) = r2_mri(1)+r2_model(1)-r2_full(1); 
   % end 
end
end 
end
disp('done.')

%% runs stats on all models 

% set rng for reproducibility 
rng(96)

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'both';

for roi = 1:size(shared_var,2)
    for layer = 1:size(shared_var,3)
        for model = 1:size(shared_var,4)

        sig_shared_var_diff(roi,layer,model) = permutation_1sample_alld(shared_var(:,roi,layer,model), nperm, cluster_th, significance_th, tail);
    end
    end 
end
[~,~,~,adj_p_shared_var_diff] = fdr_bh(sig_shared_var_diff);

%% bootstrap the peak for each model 

statsInfo.nperm = 10000;
statsInfo.cluster_th = 0.001;
statsInfo.significance_th = 0.05;
statsInfo.tail = 'both';
statsInfo.stat = [1 0]; 
nboot = 100000; 

% set rng to a fixed number 
rng(96);

for roi = 1:size(shared_var,2)
        for model = 1:size(shared_var,4)

        boot_shared_var(roi,model) = bootstrap_fixed_1D(squeeze(shared_var(:,roi,:,model)), [1:6],nboot,statsInfo);
    end
end

%% plot with error bars 


cmap = colormap('redbluetecplot');
cmap_idx = [256,200,20];
clear this_line
fig = figure;
options = [];
options.handle = fig;
options.x_axis = [1:6];
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(shared_var(:,1,:,1))*100,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(shared_var(:,1,:,2))*100,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(shared_var(:,1,:,3))*100,options);
options.color_area = cmap(ceil(20),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(20),:)%rgb('Purple');
this_line(4) = plot_areaerrorbar(squeeze(shared_var(:,1,:,4))*100,options);
for model = 1:size(shared_var,4)
    % plot stats 
    this_sig = adj_p_shared_var_diff(1,:,model);
    this_sig(this_sig>0.05) = NaN;
    this_sig(this_sig<0.05) = 1; 
    if model == 1 
        plot(options.x_axis,this_sig*-0.12*model,'.','Color','black','MarkerSize',15);
    elseif model >1
        plot(options.x_axis,this_sig*-0.12*model,'.','Color',cmap(cmap_idx(model-1),:),'MarkerSize',15);
    end 
    % plot peak CI 
    x = boot_shared_var(1,model).peak.confidence95(1);
    y = 4-0.2*model;
    err_start = 0;
    err_end = abs(boot_shared_var(1,model).peak.confidence95(2)-boot_shared_var(1,model).peak.confidence95(1));
    if model ==1, color = 'black'; elseif model >1 color = cmap(cmap_idx(model-1),:);end 
    errorbar(x, y, err_start, err_end, 'horizontal','Color',color, 'LineStyle', 'none', 'Marker', 'none', 'LineWidth', 3);

    % add observed peak value as a dot
    obs_peak_val = boot_shared_var(1,model).peak.orig;
    plot(obs_peak_val, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8)
end 
ylim([-0.5 6])
xlim([1 6])
xticks([1 3.5 6])
xticklabels({'Early';'Intermediate';'High'})
legend(this_line,'ResNet18','ResNet50','AlexNet', 'DenseNet')
title('EVC')
ylabel('Shared Variance (%)')
xlabel('Model layer')

print(fullfile(out_dir, ['dth_shared_variance_EVC.svg']), ...
             '-dsvg', '-r600')

clear this_line
fig = figure;
options = [];
options.handle = fig;
options.x_axis = [1:6];
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(shared_var(:,2,:,1))*100,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(shared_var(:,2,:,2))*100,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(shared_var(:,2,:,3))*100,options);
options.color_area = cmap(ceil(20),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(20),:)%rgb('Purple');
this_line(4) = plot_areaerrorbar(squeeze(shared_var(:,2,:,4))*100,options);
for model = 1:size(shared_var,4)
    % plot stats 
    this_sig = adj_p_shared_var_diff(2,:,model);
    this_sig(this_sig>0.05) = NaN;
    this_sig(this_sig<0.05) = 1; 
    if model == 1 
        plot(options.x_axis,this_sig*-0.12*model,'.','Color','black','MarkerSize',15);
    elseif model >1
        plot(options.x_axis,this_sig*-0.12*model,'.','Color',cmap(cmap_idx(model-1),:),'MarkerSize',15);
    end 
        % plot peak CI 
    x = boot_shared_var(2,model).peak.confidence95(1);
    y = 5.4-0.2*model;
    err_start = 0;
    err_end = abs(boot_shared_var(2,model).peak.confidence95(2)-boot_shared_var(2,model).peak.confidence95(1));
    if model ==1, color = 'black'; elseif model >1 color = cmap(cmap_idx(model-1),:);end 
    errorbar(x, y, err_start, err_end, 'horizontal','Color',color, 'LineStyle', 'none', 'Marker', 'none', 'LineWidth', 3);

    % add observed peak value as a dot
    obs_peak_val = boot_shared_var(2,model).peak.orig;
    plot(obs_peak_val, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8)
end 
ylim([-0.5 6])
xticks([1 3.5 6])
xticklabels({'Early';'Intermediate';'High'})
%legend(this_line,'ResNet18','ResNet50','AlexNet', 'DenseNet')
title('LOC')
ylabel('Shared Variance (%)')
xlabel('Model layer')

print(fullfile(out_dir, ['dth_shared_variance_LOC.svg']), ...
             '-dsvg', '-r600')

clear this_line
fig = figure;
options = [];
options.handle = fig;
options.x_axis = [1:6];
options.error = 'sem';
options.color_area = 'black';%[128 193 219]./255;    % Blue theme
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(shared_var(:,3,:,1))*100,options);
hold on
options.color_area = cmap(ceil(256),:);%rgb('DarkSeaGreen');
options.color_line = cmap(ceil(256),:);%rgb('Green');
this_line(2) = plot_areaerrorbar(squeeze(shared_var(:,3,:,2))*100,options);
options.color_area = cmap(ceil(200),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(200),:)%rgb('Purple');
this_line(3) = plot_areaerrorbar(squeeze(shared_var(:,3,:,3))*100,options);
options.color_area = cmap(ceil(20),:);%rgb('Violet');    % Orange theme
options.color_line = cmap(ceil(20),:)%rgb('Purple');
this_line(4) = plot_areaerrorbar(squeeze(shared_var(:,3,:,4))*100,options);
for model = 1:size(shared_var,4)
    % plot stats 
    this_sig = adj_p_shared_var_diff(3,:,model);
    this_sig(this_sig>0.05) = NaN;
    this_sig(this_sig<0.05) = 1; 
    if model == 1 
        plot(options.x_axis,this_sig*-0.12*model,'.','Color','black','MarkerSize',15);
    elseif model >1
        plot(options.x_axis,this_sig*-0.12*model,'.','Color',cmap(cmap_idx(model-1),:),'MarkerSize',15);
    end 
        % plot peak CI 
    x = boot_shared_var(2,model).peak.confidence95(1);
    y = 4-0.2*model;
    err_start = 0;
    err_end = abs(boot_shared_var(2,model).peak.confidence95(2)-boot_shared_var(2,model).peak.confidence95(1));
    if model ==1, color = 'black'; elseif model >1 color = cmap(cmap_idx(model-1),:);end 
    %errorbar(x, y, err_start, err_end, 'horizontal','Color',color, 'LineStyle', 'none', 'Marker', 'none', 'LineWidth', 3);

    % add observed peak value as a dot
    obs_peak_val = boot_shared_var(3,model).peak.orig;
    plot(obs_peak_val, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8)
end 
ylim([-0.5 6])
xticks([1 3.5 6])
xticklabels({'Early';'Intermediate';'High'})
legend(this_line,'ResNet18','ResNet50','AlexNet', 'DenseNet')
title('PPA')
ylabel('Shared Variance (%)')
xlabel('Model layer')


print(fullfile(out_dir, ['dth_shared_variance_PPA.svg']), ...
             '-dsvg', '-r600')


