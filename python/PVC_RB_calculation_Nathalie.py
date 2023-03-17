#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:41:36 2020

Script to calculate RBV PVC
        
    Workflow: Run function line 85 with input:
        - pet_dir = directory of pet_noPVC in seg-space
        - seg_dir = directory of segments in seg-space
        - output_pvc = directory of PVC image
        - Choose your FWHM, line 113.
        - Load directory of script in line 36.   

@author: nathaliemertens
"""

####################################################################################
#
# Packages inladen
#
####################################################################################

import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package]) 
        __import__(package)

# Load packages
packages=['os','numpy','nibabel','scipy','sys','pandas','datetime']
for pck in packages:
    import_or_install(pck)

# Load packages
import os
import numpy as np
import nibabel as nib
import scipy.ndimage.filters as filt
import scipy.optimize as opt
import sys
from tqdm import tqdm
import pandas as pd
from datetime import datetime

# Load directory of script.
#os.chdir('/Volumes/LaCie/Thomas/Projects/L3D/SCRIPTS/NATHALIE_MERTENS')

####################################################################################
#
# Pre-defined functions, used in the script
#
####################################################################################

# Function to load image data (.nii, .nii.gz)
def load_data(fname):
    nii = nib.load(fname)
#    nii = nib.as_closest_canonical(nii)
    vol = np.nan_to_num(nii.get_fdata().squeeze())
    return nii,vol,nii.affine

# Function to reshape 4D data to 3D + 1
def reshape4d(array):
    if (len(array.shape))==4:
        return array.reshape((array.shape[0]*array.shape[1]*array.shape[2],array.shape[3]))
    else:
        return array.reshape((array.shape[0]*array.shape[1]*array.shape[2],1))

####################################################################################
#
# Pre-processing
#
####################################################################################

# Define function to calculate RBV-PVC

def rb_pvc(output_prefix,pet_dir,seg_dir,output_pvc,fwhm):
            
    print('Script started.')
        
    # Define subject, used to save results
    #indx = output_pvc.find('0')
    #str_1 = output_pvc[indx-3]
    #str_2 = output_pvc[indx-2]
    #str_3 = output_pvc[indx-1]
    #str_4 = output_pvc[indx]
    #str_5 = output_pvc[indx+1]
    #str_6 = output_pvc[indx+2]
    #subject = str_1+str_2+str_3+str_4+str_5+str_6
        
    # Load data
        # PET data
    pet = load_data(pet_dir)[0]
    petData = load_data(pet_dir)[1]
    petAffine = load_data(pet_dir)[2]

    seg = load_data(seg_dir)[0]
    segData = load_data(seg_dir)[1]
    segAffine = load_data(seg_dir)[2]
        
    # Check if images are co-registered
    if np.allclose(segAffine,petAffine,rtol=1e-05):
        print('Images are co-registered.')
    else:
        print('Images are not co-registered.')
        
    ####################################################################################
    #
    # Calculate RBV PVC
    #
    ####################################################################################
        
    print('RB PVC will start.')
        
    # Optional arguments:
    #    'pet',help='PET image. Can be 4d'
    #    'seg',help='3d segmentation image'
    #    '-fwhm',help='FWHM for PET PSF in mm.
    #    '-mask',help='Mask for RBV image'
    #    '-noZero',help='Do not count 0 as a seperate ROI',action='store_const',const=1
    #    '-nii',help='Output RSF ROI averages as a ROIx1x1xTime Nifti file. Default is text file',action='store_const',const=1
    #    '-weight',help='Input previously calculated w weight matrix'
        
    #fwhm = [6.5]
    mask = None          # mask for RBV image
    noZero = 1           # 1 dus we nemen nul in de voi's niet mee
    nii = None
    weight = None
    
    #Check that images have same dimensions
    if pet.shape[0:3] != seg.shape[0:3]:
        print ('ERROR: Images do have not same dimensions...')
        sys.exit()
    
    #Setup mask
    if mask is not None:
        print('use mask')
    else:
        maskData = np.ones((seg.shape[0],seg.shape[1],seg.shape[2]))
    
    #Remove 4th dimension of segmentation
    if len(segData.shape) == 4:
        segData = segData[:,:,:,0]
    
    #Make a flattened version of the segmentation for use later
    segFlat = segData.flatten()
    
    #Reshape PET data as necessary
    if len(petData.shape) == 4:
        petData = reshape4d(petData)
        nPet = petData.shape[1]
    else:
        petData = petData.reshape((pet.shape[0]*pet.shape[1]*pet.shape[2],1))
        nPet = 1
    
    #Get ROI list and number of ROIs
    roiList = np.unique(segData)
    if noZero == 1:                         # 1, dan nemen we nul niet mee
        roiList = roiList[roiList!=0]
    nRoi = roiList.shape[0]
    
    #Make weight matrices
    if weight is not None:
        #Load w matrix
        wMatrix = np.loadtxt(weight[0])
        #Check to make sure w matrix is the correct size
        if wMatrix.shape[0] != nRoi or wMatrix.shape[1] != nRoi:
            print ('Dimensions of %s do not much the number of ROIs. Exiting...'%(weight[0]))
            sys.exit()
    else:
        wMatrix = np.zeros((nRoi,nRoi),dtype=np.float64)
    tMatrix = np.zeros((nRoi,nPet),dtype=np.float64)
    
    #Determine smoothing values for each dimension.
    voxSize = pet.header.get_zooms()
    sigmas = np.divide(fwhm[0]/np.sqrt(8.0*np.log(2.0)),voxSize[0:3])
    
    # Matrix w is de GTM matrix, berekend via:
        # = spill-over van ene voi naar andere voi
        # Step 1 = Create binary map:
            # Waarde 1 in voi i
            # Waarde 0 in alle andere voi
        # Step 2 = Calculate: binary map x psf = spill-over map
        # Step 3 = Calculate: weighted average of all other vois
    
    # Get weighted values
    for iIdx in tqdm(range(nRoi),desc='Creating SGTM Matrix'):
        # Step 1 = Maak binaire map: zet voor voi met waarde x gelijk aan 1 en elders gelijk aan 0.
        # Make i ROI. Make sure it is float as well.
        iRoi = np.float64(np.where(segData==roiList[iIdx],1.0,0.0))
        # Step 2 = Smooth i ROI for the first time, and flatten (smoothing of binary mask of roi i) = regional spread function
        iSmooth = filt.gaussian_filter(iRoi,sigmas).flatten()
        # Calculate t-values (doen voor alle tijdspunten, wij hebben er maar 1)
        # Calculate: (smoothed binary mask of roi i) x alle pet_data
        # Dit geeft dus de spill-in en spill-out weer van de pet_data door de omliggende (want binary mask is smoothed) vois.
        
        for petIdx in range(nPet):
            tMatrix[iIdx,petIdx] = iSmooth.dot(petData[:,petIdx])
        #Calculate w-matrix if necessary
        if weight is None:
            #Set diagonal weight (bereken eerst de diagonalen van w matrix)
            wMatrix[iIdx,iIdx] = iSmooth.dot(iSmooth)
            #Smooth i ROI again
            iSmooth = filt.gaussian_filter(iSmooth.reshape((seg.shape[0],seg.shape[1],seg.shape[2])),sigmas).flatten()

            #Loop thorugh other ROIs
            for jIdx in range(iIdx,nRoi):
                #Make j ROI
                jRoi = np.float64(np.where(segFlat==roiList[jIdx],1.0,0.0))
                #Get weight
                wMatrix[iIdx,jIdx] = iSmooth.dot(jRoi)
                wMatrix[jIdx,iIdx] = wMatrix[iIdx,jIdx]
    
    #Save weight matricies
    if weight is None:
        np.savetxt(os.path.join(output_pvc, output_prefix + '_wMatrix.txt'),wMatrix)
    np.savetxt(os.path.join(output_pvc,output_prefix + '_tMatrix.txt'),tMatrix)
    
    #Reshape pet data back
    if nPet==1:
        petData = petData.reshape((seg.shape[0],seg.shape[1],seg.shape[2]))
    else:
        petData = petData.reshape((seg.shape[0],seg.shape[1],seg.shape[2],nPet))
    
    #Loop through pet images
    roiCoef = np.zeros((nRoi,nPet),dtype=np.float64)
    #rbvData = np.zeros(petData.shape,dtype=np.float64)
    
    # PET_gemeten_gekend = GTM . PET_echt_nietgekend
    # b = A.x --> oplossen naar x
    # tmatrix = wmatrix . roiCoef --> oplossen naar roiCoef
    
    # First do GTM PVC --> Cgemeten = GTM x Cecht --> bereken hieruit Cecht
       
    for petIdx in tqdm(range(nPet),desc='Calculating PVC values'):
        #Get regional coefficients (PET_GTM)
        roiCoef[:,petIdx],_ = opt.nnls(wMatrix,tMatrix[:,petIdx])
        #Make Region-Spread-Function (RSF) image
        rsfData = np.zeros((seg.shape[0],seg.shape[1],seg.shape[2]),dtype=np.float64)
        # Geef elke voxel van die regio diezelfde waarde
        for roiIdx in range(nRoi):
            rsfData[segData==roiList[roiIdx]] = roiCoef[roiIdx,petIdx]

        #Make and calculate rbv image
        #rsfSmooth =  filt.gaussian_filter(rsfData,sigmas) * maskData            
        #rbvData[:,:,:,petIdx] = np.ma.divide(petData[:,:,:,petIdx] * rsfData,np.ma.array(rsfSmooth,mask=rsfSmooth==0))
    
    # (1) GTM-analysis
    #GTM-analyse = op regionaal vlak
        # => PET_av_roi = w * PET_echt_unknown
        # => PET_echt_unknown = w^-1 * PET_av_roi
    #tMatrix = PET_av_roi
    #wMatrix = w
    #roiCoef = PET_echt_unknown --> resultaat van GTM-analyse = voor elke voi hebben we een waarde 'PET_echt_unknown'
    # (2) voxel-wise analysis
    #roiCoef = s(x) = resultaat van GTM-analyse = synthetisch beeld
    #rsfData = s(x) --> Synthetisch beeld, verkregen via GTM
    #rsfSmooth = s(x)*h(x)
    #petData = f(x)
    #rbvData = f_c(x)
    
    # Save coefficients in the format user wants.
    if nii == 1:
        avg = nib.Nifti1Image(roiCoef.reshape((nRoi,1,1,nPet)),np.identity(4))
        avg.to_filename(os.path.join(output_pvc,output_prefix+ '_rsfAvg.nii'))
    else:
        np.savetxt(os.path.join(output_pvc,output_prefix + '_rsfAvg.txt'),roiCoef)
    
    # Save RBV image
    rbv_dir = os.path.join(output_pvc,output_prefix + '_PET_PVC_RB_' + str(fwhm).replace('.','_') + '_in_seg.nii')
    #rbv = nib.Nifti1Image(rbvData,petAffine)
    rbv = nib.Nifti1Image(rsfData,petAffine)
    rbv.to_filename(rbv_dir)
    print('RB PVC PET calculation done.')

    return(rbv_dir)        
    return(print('Script has finished'))
    return(rbv_dir)