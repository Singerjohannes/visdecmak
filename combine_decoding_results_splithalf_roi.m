%% combine results over iterations (for splithalf analysis) 

function combine_decoding_results_splithalf_roi(results_dir, n_perm) 

all_accs = [];

% first load the results for all iterations
for perm = 1:n_perm
    
    cfg.results.dir = fullfile(results_dir,num2str(perm));
    
    load(fullfile(cfg.results.dir, 'res_mean_decision_values.mat'))
    
    all_res_dec_val(perm,:) = results.mean_decision_values.output(:); 
        
    load(fullfile(cfg.results.dir, 'res_accuracy_minus_chance.mat'));
    
    all_accs = [results.accuracy_minus_chance.output'; all_accs]; 
    
end 

load(fullfile(cfg.results.dir, 'res_mean_decision_values.mat'))

% loop through searchlights and average results for each searchlight 
for i = 1:size(all_res_dec_val,2) 
    
    results.mean_decision_values.output{i} = mean(horzcat(all_res_dec_val{:,i}),2); 
    
end 

% save results 
save(fullfile(results_dir, 'res_mean_decision_values.mat'), 'results');

%save accuracy results
load(fullfile(cfg.results.dir, 'res_accuracy_minus_chance.mat'))

%average accuracies over perms 
results.accuracy_minus_chance.output = mean(all_accs,1); 
    
% save results 
save(fullfile(results_dir, 'res_accuracy_minus_chance.mat'), 'results');

%clear all the directories with the single results 
for perm = 1:n_perm
    
    cfg.results.dir = fullfile(results_dir,num2str(perm));
    rmdir(cfg.results.dir,'s')
    
end 
