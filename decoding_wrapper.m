%% Script for running all the decoding analyses
%
% possible analysis steps are:
% Manmade/Natural Decoding - ROI and searchlight
% Basic Level Decoding - ROI and searchlight

clear all
clc

%setup paths

parent_dir = pwd;
figure_path = fullfile(parent_dir,'figures');
betas_path = fullfile(parent_dir,'betas');
df_path = fullfile(parent_dir,'deformation_field');
behav_path = fullfile(parent_dir,'behav');
results_path = fullfile(parent_dir,'results');
roi_path = fullfile(parent_dir,'roi');
addpath(fullfile(pwd,'utils'))

% get all subject names
subs = dir(fullfile(betas_path,'*sub*'));
subs = {subs.name}';

% setup the decoding toolbox
try
    decoding_defaults;
catch
    tdt_path = input('The decoding toolbox does not seem to be on your path. Please specify the path for the folder where it is located:','s');
    addpath(tdt_path);
    % copy transres function to tdt directory, if not already there
    if ~exist(fullfile(tdt_path,'transform_results','transres_mean_decision_values.m'))
        movefile('transres_mean_decision_values.m',fullfile(tdt_path,'transform_results'));
    end
    decoding_defaults;
end

% setup spm
try
    spm;
    close all
catch
    spm_path = input('SPM does not seem to be on your path. Please specify the path of your SPM folder:','s');
    addpath(spm_path);
end


% set plot defaults

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica')

%% manmade vs. natural decoding - ROI or searchlight

for sub_idx = 1:length(subs)
    
    % select the current subject
    sub = subs{sub_idx};
    
    % set rng to fixed number for reproducibility
    rng(sub_idx);
    
    cfg = [];
    cfg.analysis = 'roi';
    cfg.n_perm = 100; %how many times should the split-half averaging and decoding be repeated
    avg_size = 2; % how many betas to average into one beta
    condition_names = cell(1,60);
    for i=1:60
        condition_names(i) = {['Image_', num2str(i)]};
    end
    beta_dir = fullfile(betas_path,sub);
    beta_avg_dir = fullfile(betas_path,sub,'avg');
    
    out_dir = fullfile(results_path,sub(1:5),'decoding','manmade_natural',cfg.analysis);
    roi_dir = fullfile(roi_path,[sub(1:5),'_rois']);
    if strcmpi(cfg.analysis, 'searchlight')
        cfg.files.mask = {fullfile(beta_dir,'mask.nii')};
    elseif strcmpi(cfg.analysis, 'roi')
        cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii');fullfile(roi_dir, 'PPA_mask.nii');fullfile(roi_dir, 'RSC_mask.nii');fullfile(roi_dir, 'OPA_mask.nii')};
    end
    
    decoding_nobetas_splithalf(condition_names,avg_size,beta_dir,beta_avg_dir,out_dir, cfg);
    
    if strcmpi(cfg.analysis,'searchlight')
        
        % combine the results from different permutations
        combine_decoding_results_splithalf(out_dir,cfg.n_perm);
        % run the distance-to-hyperplane analysis
        dth_analysis(behav_path,out_dir);
        % normalize the decoding accuracy and dth correlation maps
        
        % select files to normalize and smooth here
        fnames = cellstr(spm_select('fplistrec',fullfile(out_dir),['^[res|dth].*nii$']));
        
        spm_jobman('initcfg')
        matlabbatch = [];
        %normalize
        normparams_path = spm_select('fplist',fullfile(df_path,['sub',num2str(sub_idx,'%02i'),'_df']),['^y_struct',num2str(sub_idx,'%02i'),'.*\.(nii|img)$']); %path to forward transformation file
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {normparams_path};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fnames;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN;NaN NaN NaN];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
        
        spm_jobman('run',matlabbatch)
        
        
    elseif contains(cfg.analysis,'roi')
        combine_decoding_results_splithalf_roi(out_dir,cfg.n_perm);
        
    end
end

%% basic level decoding - ROI or searchlight

for sub_idx = 1:length(subs)
    
    % select the current subject
    sub = subs{sub_idx};
    
    % set rng to fixed number for reproducibility
    rng(sub_idx);
    
    cfg = [];
    cfg.analysis = 'roi';
    cfg.n_perm = 100; %how many times should the split-half averaging and decoding be repeated
    avg_size = 2; % how many betas to average into one beta
    condition_names = cell(1,60);
    for i=1:60
        condition_names(i) = {['Image_', num2str(i)]};
    end
    beta_dir = fullfile(betas_path,sub);
    beta_avg_dir = fullfile(betas_path,sub,'avg');
    % first average betas so that we have 2 betas for training and 1 beta for
    % testing
    n_betas = avg_betas_splithalf(condition_names,avg_size, beta_dir, beta_avg_dir,cfg);
    
    out_dir = fullfile(results_path,sub(1:5),'decoding','basic_level');
    roi_dir = fullfile(roi_path,[sub(1:5),'_rois']);
    if strcmpi(cfg.analysis, 'searchlight')
        cfg.files.mask = {fullfile(beta_dir,'mask.nii')};
    elseif strcmpi(cfg.analysis, 'roi')
        cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii');fullfile(roi_dir, 'PPA_mask.nii');fullfile(roi_dir, 'RSC_mask.nii');fullfile(roi_dir, 'OPA_mask.nii')};
    end
    
    % run the decoding for each basic-level category separately
    n_cat = 6;
    for cat_idx = 1:n_cat
        target_range = [((cat_idx - 1) * 10 + 1):cat_idx*10];
        target_condition_names = condition_names(target_range);
        if cat_idx < 4
            distractor_range = [1:30];
            distractor_condition_names = condition_names(setdiff(distractor_range,target_range));
        elseif cat_idx >3
            distractor_range = [31:60];
            distractor_condition_names = condition_names(setdiff(distractor_range,target_range));
        end
        out_dir = fullfile(results_path,sub(1:5),'decoding',['basic_level_',num2str(target_range(1)),'_',num2str(target_range(end))],cfg.analysis);
        
        decoding_nobetas_basic_level_splithalf(target_condition_names,distractor_condition_names,n_betas,beta_dir,beta_avg_dir,out_dir,cfg);
        if strcmpi(cfg.analysis,'searchlight')
            combine_decoding_results_splithalf(out_dir,cfg.n_perm);
        elseif strcmpi(cfg.analysis,'roi')
            combine_decoding_results_splithalf_roi(out_dir,cfg.n_perm);
        end 
    end
    if strcmpi(cfg.analysis,'searchlight')
        
        % combine the results from different permutations
        % combine results from different basic-level decodings and compute
        % DTH-corr
        out_dir =fullfile(results_path,sub(1:5),'decoding','basic_level');
        save_path = fullfile(results_path,sub(1:5),'decoding','basic_level','searchlight');
        combine_dth_analysis_basic_level(cfg,i_sub,out_dir,save_path);
        
        % normalize the decoding accuracy and dth correlation maps
        
        % select files to normalize and smooth here
        fnames = cellstr(spm_select('fplistrec',fullfile(out_dir),['^[res|dth].*nii$']));
        
        spm_jobman('initcfg')
        matlabbatch = [];
        %normalize
        normparams_path = spm_select('fplist',fullfile(df_path,['sub',num2str(sub_idx,'%02i'),'_df']),['^y_struct',num2str(sub_idx,'%02i'),'.*\.(nii|img)$']); %path to forward transformation file
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {normparams_path};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fnames;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN;NaN NaN NaN];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
        
        spm_jobman('run',matlabbatch)
        
    end 
end

