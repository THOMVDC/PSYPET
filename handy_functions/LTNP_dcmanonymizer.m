%% Background

% Dicom anonymizer, 
%
% Author: 
%       Thomas Vande Casteele, KU Leuven

%% work in progress, not functional yet

% Define dicom folder, dicom field and the new field value
dcmfolder='/Users/tvdcas1/Downloads/DICOM';
% field_path={'PatientName.FamilyName','PatientName.GivenName','PatientID','PatientBirthDate','OtherPatientIDs','PatientSex','PatientAge','PatientWeight'};
% field_value={'anonymized','anonymized','anonymized','anonymized','anonymized','anonymized','anonymized',0}; % has to be a number
% field_path={'PatientID','PatientBirthDate','OtherPatientIDs','PatientSex','PatientAge'};
% field_value={'anonymized','anonymized','anonymized','anonymized','anonymized'}; % has to be a number
field_path={'PatientID'};
field_value={'anonymized'}; % has to be a number

for i=1:length(field_path)
    field=field_path{i};
    value=field_value{i};
    LTNP_dcmfieldchanger(dcmfolder, field, field_value)
end

