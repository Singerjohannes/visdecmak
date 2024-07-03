%% distance to hyperplane analysis 
function combine_dth_analysis_basic_level(cfg,i_sub,res_path,out_path)

n_cat = 6;
decoding_maps =[];
for cat_idx = 1:n_cat
    target_range = [((cat_idx - 1) * 10 + 1):cat_idx*10];
    res_name = ['basic_level_',num2str(target_range(1)),'_',num2str(target_range(end))];

    load(fullfile(res_path,res_name,'searchlight','res_mean_decision_values.mat'));
    
    vol = spm_read_vols(spm_vol(fullfile(res_path,res_name,'searchlight','res_accuracy_minus_chance.nii')));
    
    decoding_maps = cat(4,decoding_maps,vol);
    all_results(cat_idx)=results;
   
end

% concatenate dec vals from different results structures 
dec_vals = cellfun(@(a,b,c,d,e,f) cat(1,a(1:10),b(1:10),c(1:10),d(1:10),e(1:10),f(1:10)),...
            all_results(1).mean_decision_values.output,...
            all_results(2).mean_decision_values.output,...
            all_results(3).mean_decision_values.output,...
            all_results(4).mean_decision_values.output,...
            all_results(5).mean_decision_values.output,...
            all_results(6).mean_decision_values.output,'UniformOutput',false);

load(fullfile('/scratch/singej96/dfg_projekt/WP1/derived', 'behav','mean_RTs_basic_categorization.mat'), 'mean_RTs')

mean_RTs = nanmean(mean_RTs);

for i = 1:length(results.mean_decision_values.output)
    
    these_dec_vals = dec_vals{i}; %reshape(results.mean_decision_values.output{i},2,60);
    if length(these_dec_vals) > 60
        these_dec_vals = mean(reshape(these_dec_vals,length(these_dec_vals)/60,60))';
    end 
    dth_corr(i) = corr(these_dec_vals,mean_RTs', 'Type','Pearson'); 
    dth_corr_only_manmade(i) = corr(these_dec_vals(1:30),mean_RTs(1:30)', 'Type','Pearson');
    dth_corr_only_natural(i) = corr(these_dec_vals(31:end),mean_RTs(31:end)', 'Type','Pearson');

end 

if contains(res_path,'normalized')
    % fill resultsvol 4D and write 4D nifi
    backgroundvalue = 0;
    % get canonical hdr from first preprocesed functional file
    template_file = dir(fullfile(cfg.sub(i_sub).dir, 'results','GLM','hrf_fitting','normalized','*.nii'));
    template_file = fullfile(template_file(1).folder,template_file(1).name);
    hdr= spm_vol(template_file); % choose canonical hdr from first classification image
    hdr = rmfield(hdr,'pinfo');
else
    %hdr = rmfield(hdr, 'dt');
    % fill resultsvol 4D and write 4D nifi
    backgroundvalue = 0;
    % get canonical hdr from first preprocesed functional file
    template_file = dir(fullfile(cfg.sub(i_sub).dir, 'alldata','run01','*.nii'));
    template_file = fullfile(template_file(1).folder,template_file(1).name);
    hdr= spm_vol(template_file); % choose canonical hdr from first classification image
    hdr = rmfield(hdr,'pinfo');
    %hdr = rmfield(hdr, 'dt');
end


if ~exist(out_path),mkdir(out_path),end 

resultsvol_hdr = hdr;
resultsvol_hdr.fname = fullfile(out_path,'res_accuracy_minus_chance.nii');
resultsvol_hdr.descrip = sprintf('Decoding accuracy averaged across tasks');
resultsvol = mean(decoding_maps,4);
spm_write_vol(resultsvol_hdr,resultsvol);
resultsvol_hdr.fname = fullfile(out_path,'dth_corr.nii');
resultsvol_hdr.descrip = sprintf('Distance to hyperplane correlation with mean subject RT map');
resultsvol = backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume with background value (default: 0)
resultsvol(results.mask_index) = dth_corr;
spm_write_vol(resultsvol_hdr,resultsvol);
resultsvol_hdr.fname = fullfile(out_path,'dth_corr_only_manmade.nii');
resultsvol(results.mask_index) = dth_corr_only_manmade;
spm_write_vol(resultsvol_hdr,resultsvol);
resultsvol_hdr.fname = fullfile(out_path,'dth_corr_only_natural.nii');
resultsvol(results.mask_index) = dth_corr_only_natural;
spm_write_vol(resultsvol_hdr,resultsvol);
end 
