% change dicom header info
dcmfolder='/Volumes/LaCie/Thomas/Projects/AGEING/UCB_MK_FLUT/RAW/UCBJ_4s_ECT/post/B079';
field_path='PatientWeight';
field_value=88; % has to be a number
%field_path='RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose';
%field_value=131090000; % has to be a number
%field_path='RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime'; % injection time
%field_path='AcquisitionTime'; %acqtime
%field_value='164037'; % Means 11:39:45 and has to be a string
LTNP_dcmfieldchanger(dcmfolder, field_path, field_value)