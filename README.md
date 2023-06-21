# visdecmak
This repository containes code for the paper "The link between visual representations and behavior in human scene perception". With the material contained in this repository all of the results in the paper can be reproduced. Link to preprint: 

## Requirements: 

To run the code in this repository you will need the following toolboxes on your matlab path: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 

The code was tested on Mac and Matlab2021a (with older Matlab versions there might be compatibility issues). 

## First-level analyses:

All first-level results can be reproduced with the decoding_wrapper.m script. 

First-level analyses are based on the beta maps for each participant. These beta maps are openly accesible via OSF (Link:). You first need to download these beta maps from OSF and then organize them with the organize_beta_maps.sh script. 

To run the decoding (ROI or searchlight) specify in the decoding_wrapper.m script which analysis you want to run (ROI or searchlight - see instructions in the decoding_wrapper.m script) and then run the script.
Depending on the type of analyses this might be time intensive (ROI results can be obtained in around 1 hour but searchlight analyses can take up to 3 to 4 hours for a single subject). 

If you select searchlight analysis, then normalization of the searchlight maps to the MNI template will be carried out. For this normalization step the deformation fields are needed for every subject. 
You need to download these files before running the decoding_wrapper.m from OSF (Link:) and then organize them with the script organize_deformation_fields.sh. 

## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.  
If you have not computed the first-level results with the steps described above you need to download the first-level results and unzip the data from OSF (Link: ) to a folder in the cloned github directory.
After you have downloaded and unzipped the files you can use group_roi_wrapper.m, group_modelling_wrapper.m and group_searchlight_wrapper.m to compute the group-level statistics for all the analyses in the paper. 
Please note that group_searchlight_wrapper.m will require big amounts of memory (up to 80GB RAM) which might preclude you from running this script on your local PC. 
In this case, you can download the precomputed group-level statistics from OSF (Link:). 