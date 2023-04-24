# visdecmak
This repository containes code for the paper "Characterizing the link between visual representations and behavior in human scene perception". With the material contained in this repository all of the results in the paper can be reproduced. Link to preprint: 

## Requirements: 

To run the code in this repository you will need the following toolboxes on your matlabpath: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 

The code was tested on Mac and Matlab2022a (with older Matlab versions there might be compatibility issues). 

## First-level analyses:

All first-level results can be reproduced with the decoding_wrapper.m script. 

First-level analyses are based on the beta maps for each participant. These beta maps are openly accesible via OSF (Link:). 
In addition, for the distance-to-hyperplane analyses you need to download the reaction time data via OSF and for normalization to the MNI template you need to download the deformation fields for each subject from OSF. 
Subsequently you need to unzip these folders and put them into a parent directory. 
To be able to run the code without changes to the paths your folder structure needs to look like the following: 

/project_dir/
/project_dir/interim_betas (folder with all the beta maps)
/poject_dir/deformation_fields 
/project_dir/behav (folder with reaction time data) 
/project_dir/visdecmak (folder with the github code stored in this repository)

If your files are organized like above you can run all the scripts from within the visdecmak folder without changing any paths in the code.  

To run the decoding (ROI or searchlight) specify in the decoding_wrapper.m script which analysis you want to run (ROI or searchlight) and then run the script.
To specify if you want to run ROI or searchlight analysis change the cfg.analysis variable to 'roi' or 'searchlight'. 
Depending on the type of analyses this might be time intensive (ROI results can be obtained in around 1 hour but searchlight analyses can take up to 3 to 4 hours for a single subject). 

If you select searchlight analysis, then normalization of the searchlight maps to the MNI template will be carried out as well. 

## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.
For the group-level analyses you need the results of the first-level analyses. You can either download them from OSF (Link:) or you can compute them yourself as described above. 

To be able to run the group-level scripts without changes to the paths your folder structure needs to look like the following: 

/project_dir/
/project_dir/results/ (folder with first-level results) 
/project_dir/visdecmak/ (folder with the github code) 

If your files are organized like above you can run all the group-level scripts from within the visdecmak folder.  

To compute the group-level results of the decoding and distance-to-hyperplane analyses for the ROIs run group_roi_wrapper.m
To compute the group-level results of the decoding and distance-to-hyperplane analyses across searchlights run group_searchlight_wrapper.m
To compute the group-level results of the commonality analysis assessing the shared variance between deep neural networks, brain and behavior run group_modelling_wrapper.m 

