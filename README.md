# PSYPET

	The psypet main pipeline takes a static PET image and a T1-weighted image as input, and produces SUV and SUVR images in native space.
	ROI-statistics from an atlas of choice are generated automatically in excel.
	The pipeline is compatible with already segmented/parcellated T1 data from other software.
	Psypet compiles handy functions such as an uptake calculator, dcm2SUV converter, centering of images on the AC-PC, ...
    Psypet is written in Matlab and relies on a correct working installation of CAT12.6.

## Main steps

	1. PET preprocessing
		- Dicom to nifti conversion
		- Centering
		- Realignment to first frame
		- Averaging all frames into one image
		- Conversion to a SUV image
	2. T1 preprocessing
		- Dicom to nifti conversion
		- Centering
		- CAT12.6 Segmentation/ Parcellation 
	3. Coregistration of processed PET to T1 image
	4. PVC if required
	5. Creation of reference region in patient space
	6. Creation of SUVR image
	7. Calculation of statistical data in ROIs from the atlas of choice


## Getting Started


### Prerequisites

	Make sure you have a correct working installation of CAT12.6
	General information: https://neuro-jena.github.io/cat/
	Download site: https://www.neuro.uni-jena.de/cat12/ -> choose cat12_r1432

### Installing

	Download the psypet code in the folder of your choice, unzip if needed, and add that folder to your Matlab search path.
	
	addpath(genpath('/folder/of/choice/PSYPET/'));

	To test this, type

	which psypet.m

	See https://www.mathworks.com/help/matlab/ref/addpath.html
	See https://www.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html

### Running 

	The main function is “psypet.m” and can be launched with “psypet_launcher.m”. 


#### Input arguments:
 
	1. subj
		- the name of the subj (string)
 
	2. T1
		- the path to the raw dicom folder, or path to the processed nifti file
        	- if raw dicom folder : we will process the dcm (conversion to nifti, centering, segmentation)
        	- if nifti file : we consider the T1 as already processed by the psypet pipeline. The atlas input argument will be considered as the segmentation and the rr as the path to the reference region in subject space (number or string)
  
	3. PET
		- the path to the raw dicom folder, or path to the processed nifti file
        	- if raw dicom folder : we will process the dcm (conversion to nifti, realignement, averaging, SUV calculation)
        	- if nifti file : we consider the PET acquisition as already processed by the psypet pipeline.
  
	4. rr
		- reference region, can be a number (5), a numeric array ([5,7,8]) or a path to mni (spm) ref region template
  
	5. atlas
        - if T1 is a nifti file, we consider "atlas" as the result of the segmentation (i.e. aparc+aseg.nii or other)
        - if T1 is a dcm folder, we consider "atlas" as the one to be obtained: 'Freesurfer', 'Fastsurfer' or 'any CAT12 suppported atlas' (ex: neuromorphometrics)

	6. outfolder
		- the output folder of your choice

#### Optional input arguments
	1. pvc (under construction)
		- 'RBV' for region based voxelwise partial volume correction
		- 'MG_orig' for Muller Gartner partial volume correction

#### Output arguments

    1. SUVR_path
		- path to the SUVR image
    2. SUVR_table_path
		- path to excel file with mean, median, uptake, volume for each ROI from the specified atlas
	3. SUV_rr_table_path
		- path to the SUV excel file of the reference region

#### Four scenarios

	1. T1 and PET data are not processed yet
		- enter the path to raw dicom files folder for both input arguments (PET, T1). Choose your CAT12 atlas preference ('neuromorphometrics','hammers',..) and psypet takes care of everything. 
	2. PET data is already processed, not T1 data
		- enter the processed nifti file path for input argument PET (psypet takes care of coregistration if not done yet) and the raw dicom files for input argument T1. Choose your atlas preference ('neuromorphometrics','hammers',..).
	3. T1 data is processed, not PET data
		- enter the processed nifti file path into the input argument T1 and the segmentation path in patient space into the input argument atlas and the reference region in patient space into the input argument rr. Raw dicom files for the input variable PET.
	4. T1 data and PET data are processed
		- see conditions of the second and third scenario's

#### General remarks

    1. Outcome images and statistics have the same voxel dimensions as the input T1

## Author

	Thomas Vande Casteele, Neuropsychiatry, KU Leuven

## Acknowledgments

	Patrick Dupont
	Michel Koole
	Nathalie Mertens
	Laura Michiels
	Maarten Laroy

## License
