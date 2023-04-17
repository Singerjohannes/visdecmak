%% distance to hyperplane analysis 
function dth_analysis(behav_path,res_path)

load(fullfile(res_path,'res_mean_decision_values.mat'));

load(fullfile(behav_path,'RT_all_subjects_5_35_categorization.mat'), 'RTs')

mean_RTs = nanmean(RTs,1); 

for i = 1:length(results.mean_decision_values.output)
    
    these_dec_vals = results.mean_decision_values.output{i}; %reshape(results.mean_decision_values.output{i},2,60);
    if length(these_dec_vals) > 60
        these_dec_vals = mean(reshape(these_dec_vals,length(these_dec_vals)/60,60))';
    end 
    dth_corr(i) = corr(these_dec_vals,mean_RTs', 'Type','Pearson'); 
    dth_corr_only_manmade(i) = corr(these_dec_vals(1:30),mean_RTs(1:30)', 'Type','Pearson');
    dth_corr_only_natural(i) = corr(these_dec_vals(31:end),mean_RTs(31:end)', 'Type','Pearson');

end 

    %hdr = rmfield(hdr, 'dt');
    % fill resultsvol 4D and write 4D nifi
    backgroundvalue = 0;
    % get canonical hdr from first preprocesed functional file
    template_file = fullfile(res_path,'res_accuracy_minus_chance.nii');
    hdr= spm_vol(template_file); % choose canonical hdr from first classification image
    hdr = rmfield(hdr,'pinfo');
    %hdr = rmfield(hdr, 'dt');

resultsvol_hdr = hdr;
resultsvol_hdr.fname = fullfile(res_path,'dth_corr.nii');
resultsvol_hdr.descrip = sprintf('Distance to hyperplane correlation with mean subject RT map');
resultsvol = backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume with background value (default: 0)
resultsvol(results.mask_index) = dth_corr;
spm_write_vol(resultsvol_hdr,resultsvol);
resultsvol_hdr.fname = fullfile(res_path,'dth_corr_only_manmade.nii');
resultsvol(results.mask_index) = dth_corr_only_manmade;
spm_write_vol(resultsvol_hdr,resultsvol);
resultsvol_hdr.fname = fullfile(res_path,'dth_corr_only_natural.nii');
resultsvol(results.mask_index) = dth_corr_only_natural;
spm_write_vol(resultsvol_hdr,resultsvol);
end 
