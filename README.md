# visdecmak
This repository containes code for the paper "Characterizing the link between visual representations and behavior in human scene perception". With the material contained in this repository all of the results in the paper can be reproduced. Link to preprint: 

## Requirements: 

To run the code in this repository you will need the following toolboxes on your matlabpath: 

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- The Decoding Toolbox 3.999 or newer (https://sites.google.com/site/tdtdecodingtoolbox/) 

The code was tested on Mac and Matlab2022a (with older Matlab versions there might be compatibility issues). 

## First-level analyses:

All first-level results can be reproduced with the decoding_wrapper.m script. 

First-level analyses are based on the beta maps for each participant. These beta maps are openly accesible via OSF (Link:). You first need to download these beta maps from OSF and then organize them with the organize_beta_maps.sh script. 

To run the decoding (ROI or searchlight) specify in the decoding_wrapper.m script which analysis you want to run (ROI or searchlight) and then run the script.
Depending on the type of analyses this might be time intensive (ROI results can be obtained in around 1 hour but searchlight analyses can take up to 3 to 4 hours for a single subject). 

If you select searchlight analysis, then normalization of the searchlight maps to the MNI template will be carried out as well. 

## Group-level analyses: 

All group-level results and the statistics in the paper can be reproduced with the code provided in this repository.  
You first need to download and unzip the data from OSF (Link: ) to a folder in the cloned github directory.
For this, download the zipped data from OSF. 
