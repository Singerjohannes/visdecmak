% function decoding_sl(condition_names,labels,beta_dir,out_dir,cfg)
%
% Wrapper script for decoding using a searchlight.
%
% Input variables:
%   condition_names: Names of all regressors to be used for classification
%   labels: Labels that should be paired with condition_names (e.g. [-1 1])
%   beta_dir: name where results are which are used for classification
%   out_dir: name of folder where results are saved
%   cfg (optional): config variable used for the decoding

function results = decoding_nobetas_splithalf(target_condition_names,distractor_condition_names,n_betas,beta_dir,beta_avg_dir,out_dir,cfg)

if ~exist('decoding_defaults.m','file')
    if ismac
        addpath('/Users/johannessinger/Documents/cloud_Berlin/Projekte/dfg/WP1/analysis_tools/tdt_3.999/decoding_toolbox');
    else
        addpath('/scratch/singej96/code_checking/reading/tdt_3.999/decoding_toolbox');
    end
end

if ~exist('cfg','var')
    cfg = decoding_defaults;
else
    cfg = decoding_defaults(cfg);
end

% Get subject number
try
    cfg.sn = str2double(regexp(beta_dir,'(?<=sub).[0-9]','match'));
catch
    warning('Could not extract subject number.')
end

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = out_dir;
cfg.results.overwrite = 1;

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
% done already

% Set the filename of your brain mask (or your ROI masks as cell matrix)
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
try cfg.files.mask;
catch
    cfg.files.mask = fullfile(beta_dir,'mask.img');
    if ~exist(cfg.files.mask,'file')
        cfg.files.mask = fullfile(beta_dir,'mask.nii');
        if ~exist(cfg.files.mask,'file')
            error('Mask not found in %s',cfg.files.mask)
        end
    end
end

% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want
% to use:

% in case similarities should be calculated
if strcmpi(cfg.decoding.software,'similarity')
    cfg.decoding.method = 'classification';
end

cfg.searchlight.unit = 'voxels'; %'mm'
cfg.searchlight.radius =4; % 12; this will yield a searchlight radius of 12mm or 10mm
cfg.searchlight.spherical = 1;
cfg.plot_design = 0;
cfg.verbose = 2; % you want all information to be printed on screen

% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q';
cfg.results.output = {'accuracy_minus_chance';'mean_decision_values'};

% Decide whether you want to see the searchlight/ROI/... during decoding
% cfg.plot_selected_voxels = 500; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% Add additional output measures if you like
% cfg.results.output = {'accuracy_minus_chance', 'AUC'}

%% Nothing needs to be changed below for a standard leave-one-run out cross
%% validation analysis.

% repeat averaging and decoding for n times
for perm = 1:cfg.n_perm
    
    % for each permutation choose a random selection of distractors to match
    % the number of targets
    rand_vec = randperm(length(distractor_condition_names),length(target_condition_names));
    these_distractor_condition_names = distractor_condition_names(rand_vec);
    condition_names = cat(2,target_condition_names,these_distractor_condition_names);
    
    % use a different folder for every iteration
    cfg.results.dir = fullfile(out_dir,num2str(perm));
    
    % Set the following field:
    % Full path to file names (1xn cell array) (e.g.
    % {'c:\exp\glm\model_button\im1.nii', 'c:\exp\glm\model_button\im2.nii', ... }
    betas = dir(fullfile(beta_avg_dir,['*',target_condition_names{1}(1:end-2),'*']));
    betas = {betas.name}';
    betas = natsortfiles(betas);
    
    % filter out only the betas corresponding to the target and distractor
    % conditions 
    for i = 1:length(betas)
        keep = false(size(condition_names));
        for j = 1:length(condition_names)
            keep(j) = strcmpi(betas{i}(1:end-15), condition_names{j});
        end
        filter_betas(i) = any(keep);
    end
    
    cfg.files.name = fullfile(beta_avg_dir,betas(filter_betas));
    % and the other two fields if you use a make_design function (e.g. make_design_cv)
    %
    if n_betas == 8
        % (1) a nx1 vector to indicate what data you want to keep together for
        % cross-validation (typically runs, so enter run numbers)
        cfg.files.chunk = repmat([1:3],1,20)';
        %
        % (2) any numbers as class labels, normally we use 1 and -1. Each file gets a
        % label number (i.e. a nx1 vector)
        cfg.files.label = [ones(1,10*3) ones(1,10*3)*2]';
        
        % This creates the leave-one-run-out cross validation design:
        cfg.design = make_design_cv(cfg);
        cfg.design.label(:,2:end) = [];
        cfg.design.set = cfg.design.set(1);
        cfg.design.train(:,1) = repmat([1 1 0],1,20);
        cfg.design.train(:,2:end) = [];
        cfg.design.test(:,1) = repmat([0 0 1],1,20);
        cfg.design.test(:,2:end)= [];
        
    elseif n_betas > 8
        % (1) a nx1 vector to indicate what data you want to keep together for
        % cross-validation (typically runs, so enter run numbers)
        cfg.files.chunk = repmat([1:4],1,20)';
        %
        % (2) any numbers as class labels, normally we use 1 and -1. Each file gets a
        % label number (i.e. a nx1 vector)
        cfg.files.label = [ones(1,10*4) ones(1,10*4)*2]';
        
        % This creates the leave-one-run-out cross validation design:
        cfg.design = make_design_cv(cfg);
        cfg.design.label(:,2:end) = [];
        cfg.design.set = cfg.design.set(1);
        cfg.design.train(:,1) = repmat([1 1 1 0],1,20);
        cfg.design.train(:,2:end) = [];
        cfg.design.test(:,1) = repmat([0 0 0 1],1,20);
        cfg.design.test(:,2:end)= [];
    end
    
    % run decoding
    results = decoding(cfg);
end
end