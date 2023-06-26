# visdecmak
This repository containes code for the paper "The link between visual representations and behavior in human scene perception". With the material contained in this repository all of the results in the paper can be reproduced. Link to preprint: 

## Requirements: 

To run the code in this repository you will need the following toolboxes on your matlab path: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 

The code was tested on Mac and Matlab2021a (with older Matlab versions there might be compatibility issues). 

## First-level analyses:

All first-level results can be reproduced with the decoding_wrapper.m script. 

First-level analyses are based on the beta maps for each participant. These beta maps, the subject-specific ROI masks and the subject-specific deformation fields (for normalizing the results to the MNI template) are openly accesible via OSF (Link:). You first need to download the fMRI data component from OSF and then unzip and move them with the organize_fMRI_data.sh script. 

To run the decoding (ROI or searchlight) specify in the decoding_wrapper.m script which analysis you want to run (ROI or searchlight - see instructions in the decoding_wrapper.m script) and then run the script.
Depending on the type of analyses this might be time intensive (ROI results can be obtained in around 1 hour but searchlight analyses can take up to 3 to 4 hours for a single subject). 

## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.  

If you have not computed the first-level results with the steps described above you need to download the data from the first-level results component on OSF and unzip and move the data with the organize_first_level_results.sh script. 

After you have downloaded and unzipped the files you can use group_roi_wrapper.m, group_modelling_wrapper.m and group_searchlight_wrapper.m to compute the group-level statistics for all the analyses in the paper. 

Please note that group_searchlight_wrapper.m will require big amounts of memory (up to 80GB RAM) which might preclude you from running this script on your local PC. 
In this case, you can download the precomputed group-level statistics from the group-level results component on OSF. 
