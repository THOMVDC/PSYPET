%% Background

% Launcher example for LTNP_dcmfieldchanger
%
% Author: 
%       Thomas Vande Casteele, KU Leuven

%% Settings

% Define dicom folder, dicom field and the new field value
dcmfolder='/Volumes/LaCie/Thomas/Projects/RETINAL_IMAGING/DATA/RAW/AMYLOID/MCI009';
field_path='PatientWeight';
field_value=88; % has to be a number

% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose';
% field_value=131090000; % has to be a number
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime'; % injection time
% field_path='AcquisitionTime'; %acqtime
% field_value='164037'; % Means 11:39:45 and has to be a string
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose';
% field_value=125860000; % has to be a number
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime';
% field_value='113945' % Means 11:39:45 and has to be a string
% field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife';
% field_value=1.2231e+03
% field_path='AcquisitionTime';
% field_value='164037'; % Means 16:40:37 and has to be a string

%% Launch for single subject
LTNP_dcmfieldchanger(dcmfolder, field_path, field_value)