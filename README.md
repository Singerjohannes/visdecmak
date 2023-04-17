# visdecmak
This repository containes code for the paper "Characterizing the link between visual representations and behavior in human scene perception". With the material contained in this repository all of the results in the paper can be reproduced. Link to preprint: 

## Requirements: 

To run the code in this repository you will need the following toolboxes on your matlabpath: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 

The code was tested on Mac and Matlab2022a (with older Matlab versions there might be compatibility issues). 

## First-level analyses:

We provide beta maps for every subject to demonstrate how the first-level results are computed.  

To run the decoding (ROI or searchlight) run decoding_wrapper.m and specify in the script which analysis you want to run (ROI or searchlight).
Depending on the type of analyses this might be time intensive (ROI results can be obtained in around 1 hour but searchlight analyses can take up to 3 to 4 hours for a single subject). 

## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.  
You first need to download and unzip the data from OSF (Link: https://osf.io/vsc6y/) to a folder in the cloned github directory.
For this, download the zipped data from the MEG data and fMRI data components on OSF separately. 
Then run the scripts "unzip_move_fmri_data.sh" and "unzip_move_meg_data.sh" to unzip and move the data to the following folder structure: 