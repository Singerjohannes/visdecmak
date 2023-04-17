%% Script for running all the decoding analyses
%
% possible analysis steps are:
% Manmade/Natural Decoding - ROI and searchlight

clear all
clc

%setup paths

path = pwd;
figure_path = fullfile(path,'figures');
betas_path = '/Users/johannessinger/Documents/cloud_Berlin/Projekte/dfg/WP1/paper/data/betas';
behav_path = '/Users/johannessinger/Documents/cloud_Berlin/Projekte/dfg/WP1/paper/data';
results_path = '/Users/johannessinger/Documents/cloud_Berlin/Projekte/dfg/WP1/paper/results';
roi_path = '/Users/johannessinger/Documents/cloud_Berlin/Projekte/dfg/WP1/paper/data/roi';

% get all subject names
subs = dir(fullfile(betas_path,'*sub*'));
subs = {subs.name}';

% add utils

addpath(fullfile(path,'utils'));

% add first level functions

addpath(fullfile(path,'first_level','fmri'));

% setup the decoding toolbox
try
    decoding_defaults;
catch
    tdt_path = input('The Decoding Toolbox seems to be not on your path. Please enter the path to your TDT version:\n','s');
    addpath(tdt_path);
    decoding_defaults;
end

% setup spm
try
    spm;
    close all
catch
    spm_path = input('SPM seems to be not on your path. Please enter the path to your SPM version:\n','s');
    addpath(spm_path);
end


% set plot defaults

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3, 'DefaultAxesFontName', 'Helvetica','DefaultTextFontName', 'Helvetica')

% get colormap
%cmap = colormap('redblueTecplot');
close all

%% manmade vs. natural decoding - ROI or searchlight

for sub_idx = 1:length(subs)
    
    % select the current subject
    sub = subs{sub_idx};
    
    cfg = [];
    cfg.analysis = 'roi';
    cfg.n_perm = 2; %how many times should the split-half averaging and decoding be repeated
    avg_size = 2 ; % how many betas to average into one beta
    condition_names = cell(1,60);
    for i=1:60
        condition_names(i) = {['Image_', num2str(i)]};
    end
    beta_dir = fullfile(betas_path,sub);
    beta_avg_dir = fullfile(betas_path,sub,'avg');
    
    out_dir = fullfile(results_path,sub(end-4:end),'decoding','manmade_natural',cfg.analysis);
    roi_dir = fullfile(roi_path,sub(end-4:end));
    if strcmpi(cfg.analysis, 'searchlight')
        cfg.files.mask = {fullfile(beta_dir,'mask.nii')};
    elseif strcmpi(cfg.analysis, 'roi')
        cfg.files.mask = {fullfile(roi_dir, 'evcmask.nii');fullfile(roi_dir, 'loc_mask.nii');fullfile(roi_dir, 'PPA_mask.nii')};
    end
    
    decoding_nobetas_splithalf(condition_names,avg_size,beta_dir,beta_avg_dir,out_dir, cfg);
    
    if strcmpi(cfg.analysis,'searchlight')
        
        combine_decoding_results_splithalf(out_dir,cfg.n_perm);
        dth_analysis(behav_path,out_dir);
        
    elseif contains(cfg.analysis,'roi')
        combine_decoding_results_splithalf_roi(out_dir,cfg.n_perm);
        
    end
end
