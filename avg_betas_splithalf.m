%% this function averages the beta images for later RSA 
function n_betas = avg_betas_splithalf(labelnames, avg_size,beta_dir,out_dir,cfg)

% check if avg dir already exists 
if exist(out_dir)
    delete([out_dir, '/*.nii'])
end 

if ~exist(out_dir), mkdir(out_dir), end; 

% since the labels are arbitrary, we will set them randomly to -1 and 1
labels = [ones(1,30) ones(1,30)*2];

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% now retrieve number of runs from regressor names 
n_betas = max(cell2mat(regressor_names(2,:)));

% Extract all information for the cfg.files structure (labels will be [1 -1] )
cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir);

% create the shuffle vector before so its the same for all conditions 
shuffle_vector = randperm(n_betas);

for this_cat = 1: length(labelnames) 
    
    these_files = cfg.files.name(find(ismember(cfg.files.labelname, labelnames{this_cat})));
    
    % shuffle the files 
    these_files = these_files(shuffle_vector);
    
    beta_all = [];
    
    % load all betas in one matrix 
    for this_file = 1:length(these_files) 
        
    vol = spm_vol(these_files{this_file});
    beta = spm_read_vols(vol);
    
    beta_all = cat(4,beta_all, beta); 
    
    end 
    
    %average 2 betas into one beta for the first half 
    ct = 1; 
    for i=1:avg_size:length(these_files)/2 
        
        beta_avg = mean(beta_all(:,:,:,i:i+(avg_size-1)),4); 
    
    vol.fname = fullfile(out_dir,[labelnames{this_cat},'_beta_avg_' num2str(ct) '.nii']);
    spm_write_vol(vol,beta_avg);
    ct=ct+1; 
    end 
    ct = ct+1; 
    % average all betas into one for the second half 
    beta_avg = mean(beta_all(:,:,:,(i+2):end),4); 
    
    vol.fname = fullfile(out_dir,[labelnames{this_cat},'_beta_avg_' num2str(ct) '.nii']);
    spm_write_vol(vol,beta_avg);
    
end        
end 