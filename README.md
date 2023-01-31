# PSYPET

The psypet processing pipeline is written in Matlab and relies on a correct working installation of CAT12.6 (https://neuro-jena.github.io/cat/).

Main steps:
	1/ PET preprocessing
		Dicom to nifti conversion
		Centering
		Avering to first frame
		Realignment to first frame
		Avering of all frames into one image
		Conversion to a SUV image
	2/ T1 processing
		Dicom to nifti conversion
		Centering
		Segmentation
		Parcellation 
	3/ Coregistration of processed PET and T1 images
	4/ PVC RBV
	5/ Creation of reference region in patient space
	6/ Creation of SUVR image
	7/ Calculation of statistical data in ROIs


The pipeline is compatible with already segmented/parcellated T1 data from other software.


## Getting Started


### Prerequisites

	Make sure you have a correct working installation of CAT12.6

### Installing

	Download the software in the folder of your choice, unzip if needed, and add the folder to your Matlab search path.
	See https://www.mathworks.com/help/matlab/ref/addpath.html
	See https://www.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html

## Running 

The main function is “psypet.m” and can be launched with “psypet_launcher.m”. It takes following arguments as input:
 
	1/ subj : the name of the subj (string)
 
	2/ T1 : the path to the raw dicom folder, or path to the processed nifti file
               * if raw dicom folder : we will process the dcm (conversion
               to nifti, centering, segmentation)
               * if nifti file : we consider the T1 as already processed
               by the psypet pipeline. The atlas input argument will be
               considered as the segmentation and the rr as the path to
               the reference region in subject space (number or string)
  
	3/ PET :  the path to the raw dicom folder, or path to the processed nifti file
               * if raw dicom folder : we will process the dcm (conversion
               to nifti, realignement, averaging, SUV calculation)
               * if nifti file : we consider the PET acquisition as already processed
               by the psypet pipeline.
  
	4/ rr : reference region, can be a number (5), a numeric array ([5,7,8]) or a path to mni (spm) ref region
   template
  
	5/ atlas/segmentation :
               * if T1 is a nifti file, we consider "atlas" as the result of the segmentation (i.e. aparc+aseg.nii or other)
                * if T1 is a dcm folder, we consider "atlas" as the one to be obtained: 'Freesurfer', 'Fastsurfer' or 'any CAT12 suppported atlas'
   (ex: neuromorphometrics)

	6/ outfolder : the output folder of your choice

## Output

The expected output:

	1/ SUVR image, SUVR PVC image
	2/ Excel file with mean, median, uptake, volume for each ROI from the specified atlas


## Authors

## License

## Acknowledgments

