%% group searchlight wrapper
% this script contains all the group analyses for the searchlight results 
% 1. Manmade/Natural or Basic-level Decoding + Statistics 
% 2. Distance-to-hyperplane correlation with the categorization RTs +
% statistics 
% 3. Distance-to-hyperplane correlation with the distraction RTs +
% statistics

clear all 
clc

%setup paths 

path = pwd;
figure_path = fullfile(path,'figures');
% create figure path
if ~isdir(figure_path); mkdir(figure_path); end 
% specify path where first level results are stored
results_path = fullfile(path, 'results'); 
% specify where group level results should be stored
out_dir = fullfile(results_path,'group');

% add stats functions 

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

%% decoding + statistics

% get fmri subnames 

fmri_subs = dir(fullfile(results_path,'sub*'));
fmri_subs = {fmri_subs.name}';

% specify excluded subjects

fmri_excluded_subs = {'sub12'};  %sub12 was excluded because didnt fit inclusion criteria

decoding_maps = [];

%specify results name
res_name = 'manmade_natural'; % change to "basic-level" in case you want to compute statistics for basic-level decoding 
fname = 'wres_accuracy_minus_chance.nii';

% load searchlight results

for sub_no = 1:length(fmri_subs)
    
    sub_id = fmri_subs{sub_no};
    
    if ~any(ismember(fmri_excluded_subs,sub_id))
        
        % load fMRI RDMs
        fmri_fname = fullfile(results_path,sub_id,'decoding',res_name,'searchlight',fname);

        if exist(fmri_fname)
            fprintf('Loading fMRI %s\n',sub_id);
            
            vol =spm_read_vols(spm_vol(fmri_fname));
            vol(isnan(vol))= 0; % convert NaNs to 0s
            decoding_maps = cat(4, decoding_maps, vol);

        else
            fprintf('Results not complete for sub %s\n',sub_id);
            
        end
        
    end
end

% find for each voxel how many subjects have values 

sz = size(decoding_maps);
stats_mask = zeros(sz(1),sz(2),sz(3));

for vox = 1:(sz(1)*sz(2)*sz(3))
    
    
   [idx1,idx2,idx3] = ind2sub([sz(1),sz(2),sz(3)],vox);
   
   stats_mask(vox) = sum(decoding_maps(idx1,idx2,idx3,:) ~=0);
   
end

% take only those voxels for which more than 70% of the subjects have
% values
cutoff = ceil(size(decoding_maps,4)*0.7);
stats_mask = stats_mask>cutoff; 

% now mask the group results with the stats mask for running the stats 

decoding_maps_masked = [];
for sub= 1:size(decoding_maps,4)
    
      this_vol = decoding_maps(:,:,:,sub);
      this_vol(stats_mask==0) = NaN;
      decoding_maps_masked = cat(4,decoding_maps_masked,this_vol);
end 

% write mean results without stats masking
hdr = spm_vol(fmri_fname); %  this should be an image with the same dimensionality as the searchlight results
%hdr = rmfield(hdr, 'dt'); % get rid of scaling factors from the original image
hdr.descrip = sprintf('Searchlight decoding');
%hdr = rmfield(hdr, 'n');
hdr.fname = fullfile(out_dir,[res_name,'_accuracy_minus_chance.nii']); %accuracy_minus_chance
mean_vol = mean(decoding_maps,4);
spm_write_vol(hdr, mean_vol);
      
% compute stats 

%set rng to fixed number for reproducibility 
rng(96); 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'right';

[sig_searchlight_max,sig_searchlight_weigh,~,~,~,statmap] = permutation_cluster_1sample_weight_alld_less_mem (permute(decoding_maps_masked,[4 1 2 3]), nperm, cluster_th, significance_th, tail);

% write results 
hdr = spm_vol(fmri_fname); %  this should be an image with the same dimensionality as the searchlight results
hdr = rmfield(hdr, 'dt'); % get rid of scaling factors from the original image
hdr.descrip = sprintf('Searchlight decoding masked');
hdr = rmfield(hdr, 'n');
hdr.fname = fullfile(out_dir,[res_name,'_accuracy_minus_chance_masked_max.nii']);
mean_vol = mean(decoding_maps,4);
spm_write_vol(hdr, mean_vol.*sig_searchlight_max);


%% dth correlation with categorization RTs + statistics 

decoding_maps = [];

%specify results name
res_name = 'manmade_natural'; %change to 'basic-level' in case you want to run stats on the basic-level results
fname = 'wdth_corr.nii'; %'s05wres_accuracy_minus_chance.nii' ;

% load searchlight results

for sub_no = 1:length(fmri_subs)
    
    sub_id = fmri_subs{sub_no};
    
    if ~any(ismember(fmri_excluded_subs,sub_id))
        
        % load fMRI RDMs
        fmri_fname = fullfile(results_path,sub_id,'decoding',res_name,'searchlight',fname);
        
        if exist(fmri_fname)
            fprintf('Loading fMRI %s\n',sub_id);
            
            hdr = spm_vol(fmri_fname);
            vol = spm_read_vols(hdr );
            decoding_maps = cat(4, decoding_maps, vol);
        else
            fprintf('Results not complete for sub %s\n',sub_id);
            
        end
        
    end
end

% now mask the group results with the stats mask for running the stats 
decoding_maps_masked = [];
for sub= 1:size(decoding_maps,4)
    
      this_vol = decoding_maps(:,:,:,sub);
      this_vol(stats_mask==0) = NaN;
      decoding_maps_masked = cat(4,decoding_maps_masked,this_vol);
end 

% compute stats 

%set rng to fixed number for reproducibility 
rng(96); 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'both';

[sig_searchlight_max,sig_searchlight_weigh,~,~,~,statmap] = permutation_cluster_1sample_weight_alld_less_mem (permute(decoding_maps_masked,[4 1 2 3]), nperm, cluster_th, significance_th, tail);

hdr = spm_vol(fmri_fname); %  this should be an image with the same dimensionality as the searchlight results
hdr = rmfield(hdr, 'dt'); % get rid of scaling factors from the original image
hdr.descrip = sprintf('Distance-to-hyperplane correlation masked');
hdr = rmfield(hdr, 'n');

hdr.fname = fullfile(out_dir,[res_name,'_dth_corr_masked_max.nii']);
mean_vol = mean(decoding_maps,4);
spm_write_vol(hdr, mean_vol.*sig_searchlight_max);

%% dth correlation with distraction RTs + statistics 

decoding_maps = [];

%specify results name
res_name = 'manmade_natural';
fname = 'wdth_corr_distraction.nii'; %'s05wres_accuracy_minus_chance.nii' ;

% load searchlight results

for sub_no = 1:length(fmri_subs)
    
    sub_id = fmri_subs{sub_no};
    
    if ~any(ismember(fmri_excluded_subs,sub_id))
        
        % load fMRI RDMs
        fmri_fname = fullfile(results_path,sub_id,'decoding',res_name,'searchlight',fname);
        
        if exist(fmri_fname)
            fprintf('Loading fMRI %s\n',sub_id);
            
            hdr = spm_vol(fmri_fname);
            vol = spm_read_vols(hdr );
            decoding_maps = cat(4, decoding_maps, vol);
        else
            fprintf('Results not complete for sub %s\n',sub_id);
            
        end
        
    end
end

% now mask the group results with the stats mask for running the stats 
decoding_maps_masked = [];
for sub= 1:size(decoding_maps,4)
    
      this_vol = decoding_maps(:,:,:,sub);
      this_vol(stats_mask==0) = NaN;
      decoding_maps_masked = cat(4,decoding_maps_masked,this_vol);
end 

% compute stats 

%set rng to fixed number for reproducibility 
rng(96); 

nperm = 10000;
cluster_th = 0.001;
significance_th = 0.05;
tail = 'both';

[sig_searchlight_max,sig_searchlight_weigh,~,~,~,statmap] = permutation_cluster_1sample_weight_alld_less_mem (permute(decoding_maps_masked,[4 1 2 3]), nperm, cluster_th, significance_th, tail);

hdr = spm_vol(fmri_fname); %  this should be an image with the same dimensionality as the searchlight results
hdr = rmfield(hdr, 'dt'); % get rid of scaling factors from the original image
hdr.descrip = sprintf('Distance-to-hyperplane correlation masked');
hdr = rmfield(hdr, 'n');
hdr.fname = fullfile(out_dir,[res_name,'_dth_corr_distraction_masked_max.nii']);
mean_vol = mean(decoding_maps,4);
spm_write_vol(hdr, mean_vol.*sig_searchlight_max);


