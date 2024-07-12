%% group modelling wrapper
% this script contains the group analyses for the neural network modelling 
% for basic level categorization

clear all 
clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');
% create figure path
if ~isdir(figure_path); mkdir(figure_path); end 
% specify path where first level results are stored
results_path = fullfile(path, 'results'); 
% specify path where neural network modelling results are stored 
modelling_path = fullfile(path,'modelling','basic_level');
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
    spm_path = input('Please insert the path to your SPM folder:');
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

%% load fmri and behavioral results for basic-level categorization

% get fmri subnames 

fmri_subs = dir(fullfile(results_path,'*sub*'));
fmri_subs = {fmri_subs.name}';
fmri_subs = fmri_subs(1:end-1);

% specify some constants
n_cat = 6;
n_roi = 3;

% specify excluded subjects
excluded_subjects = {'sub12'};
% intialiaze results
all_dec_vals = cell(length(fmri_subs)-1,n_roi);
all_decoding_roi = [];

for cat_idx = 1:n_cat
    
    target_range = [((cat_idx - 1) * 10 + 1):cat_idx*10];
    
    %specify results name
    res_name = ['basic_level_',num2str(target_range(1)),'_',num2str(target_range(end))];
    disp(res_name)
    decoding_roi = [];
    
    for i_sub = 1:length(fmri_subs)-1
        
        sub_id = fmri_subs{i_sub};
        
        if any(ismember(excluded_subjects, sub_id)), continue, end
        
        if ~isdir(fullfile(results_path,sub_id,'decoding',res_name,'roi')); fprintf('Not %i\n',i_sub);  continue; end;
        
        results_dir =  fullfile(results_path,sub_id,'decoding',res_name,'roi');
        load(fullfile(results_dir, 'res_accuracy_minus_chance.mat'));
        decoding_roi = [results.accuracy_minus_chance.output; decoding_roi];
        load(fullfile(results_dir, 'res_mean_decision_values.mat'));
        for roi = 1:size(decoding_roi,2)
            these_dec_vals = results.mean_decision_values.output{roi};
            all_dec_vals(i_sub,roi) = {[all_dec_vals{i_sub,roi};these_dec_vals(1:10)]};
        end
    end
    
    fprintf('Mean decoding accuracies over all subjects EVC: %2f, LOC: %2f, PPA: %2f\n', mean(decoding_roi(:,1)),mean(decoding_roi(:,2)),mean(decoding_roi(:,3)));
    all_decoding_roi = cat(3,all_decoding_roi,decoding_roi);
   
end

% average decoding across basic-level categories
decoding_roi = mean(all_decoding_roi,3);

%exclude subject 12 from decision value cell array
all_dec_vals(12,:) = [];

% convert decision values to matrix

dec_vals_mat = [];
for i_sub = 1:size(all_dec_vals,1)
    
    dec_vals_mat(i_sub,:,:) = cell2mat(all_dec_vals(i_sub,:));
end



% load behavior
load(fullfile(behav_path,'mean_RTs_basic_categorization.mat'))

mean_RTs = nanmean(mean_RTs,1);%.*[ones(1,30), ones(1,30)*-1]


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

order_idx = [4 5 7 8 6 3]; % ordered values from low to high level features
alexnet_distances =alexnet_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

% load densenet distances

dist_files = dir(fullfile(modelling_path,'*.mat'));
dist_filenames = {dist_files.name}'; 
filenames_info = regexpi(dist_filenames,'^densenet161_places.*'); 

densenet_distances = [];
for idx = 1:length(dist_filenames)
    
    if filenames_info{idx}
    results = load(fullfile(modelling_path,dist_filenames{idx}));
    densenet_distances = cat(2,densenet_distances,abs(results.dec_vals'));
    end 
end 

order_idx = [6 1 2 3 4 5]; % ordered values from low to high level features
densenet_distances =densenet_distances(:,order_idx);

% little check that the order idxs are correct 
check_filenames = dist_filenames(find(~cellfun(@isempty,filenames_info)));
check_filenames = check_filenames(order_idx);
disp(check_filenames)

%% shared variance for each model independently - ResNet18

disp('Running commonality analysis...')

shared_var = [];
shared_var_model = [];
shared_var_mri = [];
for sub = 1:size(dec_vals_mat,1)
for j_roi = 1:3
    for model_idx = 1:size(resnet18_distances,2)
        
    xMRI = squeeze(dec_vals_mat(sub,:,j_roi))'; 
    xmodel = resnet18_distances(:,model_idx);
    
    y = mean_RTs';
    
    % calculate regressions
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);
    
    % substract according to formula
    shared_var_mri(sub,j_roi) = r2_mri(1);
    shared_var(sub,j_roi,model_idx,1) = r2_mri(1)+r2_model(1)-r2_full(1);
end
end 
end
disp('done.')

%% shared variance for each model independently - ResNet50

disp('Running commonality analysis...')


for sub = 1:size(dec_vals_mat,1)
for j_roi = 1:3
    for model_idx = 1:size(resnet50_distances,2)
        
    xMRI = squeeze(dec_vals_mat(sub,:,j_roi))'; 
    xmodel = resnet50_distances(:,model_idx);
    
    y = mean_RTs';
    
    % calculate regressions
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);
    
    % substract according to formula
    shared_var(sub,j_roi,model_idx,2) = r2_mri(1)+r2_model(1)-r2_full(1);

end
end 
end
disp('done.')

%% shared variance for each model independently - AlexNet

disp('Running commonality analysis...')

for sub = 1:size(dec_vals_mat,1)
for j_roi = 1:3
    for model_idx = 1:size(alexnet_distances,2)
        
    xMRI = squeeze(dec_vals_mat(sub,:,j_roi))'; 
    xmodel = alexnet_distances(:,model_idx);
    
    y = mean_RTs';
    
    % calculate regressions
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    % substract according to formula
    shared_var(sub,j_roi,model_idx,3) = r2_mri(1)+r2_model(1)-r2_full(1); 
end
end 
end
disp('done.')

%% shared variance for each model independently - DenseNet

disp('Running commonality analysis...')

for sub = 1:size(dec_vals_mat,1)
for j_roi = 1:3
    for model_idx = 1:size(densenet_distances,2)
        
    xMRI = squeeze(dec_vals_mat(sub,:,j_roi))'; 
    xmodel = densenet_distances(:,model_idx);
    
    y = mean_RTs';
    
    % calculate regressions
    [~,~,~,~,r2_mri] = regress(y, [ones(1,60)' xMRI]);
    [~,~,~,~,r2_model] = regress(y,[ones(1,60)' xmodel]);
    [~,~,~,~,r2_full] = regress(y, [ones(1,60)' xmodel xMRI]);

    % substract according to formula
    shared_var(sub,j_roi,model_idx,4) = r2_mri(1)+r2_model(1)-r2_full(1); 

    end
end 
end
disp('done.')

%% runs stats on all models 

% set rng for reproducibility 
rng(96,'twister')

% set stats defaults 
nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

for roi = 1:2 % calculate stats only for EVC and LOC
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
statsInfo.tail = 'right';
statsInfo.stat = [1 0]; 
nboot = 100000; 

% set rng to a fixed number 
rng(96,'twister');

for roi = 1:2
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
options.color_area = cmap(ceil(256),:);
options.color_line = cmap(ceil(256),:);
this_line(2) = plot_areaerrorbar(squeeze(shared_var(:,1,:,2))*100,options);
options.color_area = cmap(ceil(200),:);
options.color_line = cmap(ceil(200),:);
this_line(3) = plot_areaerrorbar(squeeze(shared_var(:,1,:,3))*100,options);
options.color_area = cmap(ceil(20),:);
options.color_line = cmap(ceil(20),:);
this_line(4) = plot_areaerrorbar(squeeze(shared_var(:,1,:,4))*100,options);
for model = 1:size(shared_var,4)
    % plot stats 
    this_sig = adj_p_shared_var_diff(1,:,model);
    this_sig(this_sig>0.05) = NaN;
    this_sig(this_sig<0.05) = 1; 
    if model == 1 
        plot(options.x_axis,this_sig*-0.12*model-0.5,'.','Color','black','MarkerSize',15);
    elseif model >1
        plot(options.x_axis,this_sig*-0.12*model-0.5,'.','Color',cmap(cmap_idx(model-1),:),'MarkerSize',15);
    end 
    % plot peak CI 
    x = boot_shared_var(1,model).peak.confidence95(1);
    y = 4.5-0.2*model;
    err_start = 0;
    err_end = abs(boot_shared_var(1,model).peak.confidence95(2)-boot_shared_var(1,model).peak.confidence95(1));
    if model ==1, color = 'black'; elseif model >1 color = cmap(cmap_idx(model-1),:);end 
    errorbar(x, y, err_start, err_end, 'horizontal','Color',color, 'LineStyle', 'none', 'Marker', 'none', 'LineWidth', 3);

    % add observed peak value as a dot
    obs_peak_val = boot_shared_var(1,model).peak.orig;
    plot(obs_peak_val, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8)
end 
x = 1:6;
y1 = ones(1,6)*(mean(shared_var_mri(:,1)) - std(shared_var_mri(:,1)) / sqrt(size(shared_var_mri,1))) * 100;
y2 = ones(1,6)*(mean(shared_var_mri(:,1)) + std(shared_var_mri(:,1)) / sqrt(size(shared_var_mri,1))) * 100;
% Filling the area between y1 and y2
patch([x, fliplr(x)], [y1, fliplr(y2)],[0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
this_line(end+1) = plot(x,ones(1,6)*mean(shared_var_mri(:,1))*100,'Color',[0.5 0.5 0.5])
ylim([-1 4.5])
xlim([1 6])
xticks([1 3.5 6])
xticklabels({'Early';'Intermediate';'High'})
%legend(this_line,'ResNet18','ResNet50','AlexNet', 'DenseNet','Brain-behavior')
title('EVC')
ylabel('Shared Variance (%)')
xlabel('Model layer')

print(fullfile(out_dir, ['dth_shared_variance_basic_level_EVC.svg']), ...
             '-dsvg', '-r600')

%% plot for LOC
fig = figure;
clear this_line
options = [];
options.handle = fig;
options.x_axis = [1:6];
options.error = 'sem';
options.color_area = 'black';
options.color_line = [17 17 17]./255;
options.alpha      = 0.5;
options.line_width = 3;
this_line(1) = plot_areaerrorbar(squeeze(shared_var(:,2,:,1))*100,options);
hold on
options.color_area = cmap(ceil(256),:);
options.color_line = cmap(ceil(256),:);
this_line(2) = plot_areaerrorbar(squeeze(shared_var(:,2,:,2))*100,options);
options.color_area = cmap(ceil(200),:);
options.color_line = cmap(ceil(200),:);
this_line(3) = plot_areaerrorbar(squeeze(shared_var(:,2,:,3))*100,options);
options.color_area = cmap(ceil(20),:);
options.color_line = cmap(ceil(20),:);
this_line(4) = plot_areaerrorbar(squeeze(shared_var(:,2,:,4))*100,options);
for model = 1:size(shared_var,4)
    % plot stats 
    this_sig = adj_p_shared_var_diff(2,:,model);
    this_sig(this_sig>0.05) = NaN;
    this_sig(this_sig<0.05) = 1; 
    if model == 1 
        plot(options.x_axis,this_sig*-0.12*model-0.5,'.','Color','black','MarkerSize',15);
    elseif model >1
        plot(options.x_axis,this_sig*-0.12*model-0.5,'.','Color',cmap(cmap_idx(model-1),:),'MarkerSize',15);
    end 
    % plot peak CI 
    x = boot_shared_var(1,model).peak.confidence95(1);
    y = 4-0.2*model;
    err_start = 0;
    err_end = abs(boot_shared_var(2,model).peak.confidence95(2)-boot_shared_var(2,model).peak.confidence95(1));
    if model ==1, color = 'black'; elseif model >1 color = cmap(cmap_idx(model-1),:);end 
    errorbar(x, y, err_start, err_end, 'horizontal','Color',color, 'LineStyle', 'none', 'Marker', 'none', 'LineWidth', 3);

    % add observed peak value as a dot
    obs_peak_val = boot_shared_var(2,model).peak.orig;
    plot(obs_peak_val, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerSize', 8)
end 
x = 1:6;
y1 = ones(1,6)*(mean(shared_var_mri(:,2)) - std(shared_var_mri(:,2)) / sqrt(size(shared_var_mri,1))) * 100;
y2 = ones(1,6)*(mean(shared_var_mri(:,2)) + std(shared_var_mri(:,2)) / sqrt(size(shared_var_mri,1))) * 100;
% Filling the area between y1 and y2
patch([x, fliplr(x)], [y1, fliplr(y2)],[0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
this_line(end+1) = plot(x,ones(1,6)*mean(shared_var_mri(:,2))*100,'Color',[0.5 0.5 0.5])
ylim([-1  4.5])
xlim([1 6])
xticks([1 3.5 6])
xticklabels({'Early';'Intermediate';'High'})
%legend(this_line,'ResNet18','ResNet50','AlexNet', 'DenseNet','Brain-behavior')
title('LOC')
ylabel('Shared Variance (%)')
xlabel('Model layer')

print(fullfile(out_dir, ['dth_shared_variance_basic_level_LOC.svg']), ...
             '-dsvg', '-r600')