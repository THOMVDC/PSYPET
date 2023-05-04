function [ptname,ptID,bodyweight, injdosis, acqtime, injtime, halftime,rescale_factor]=LTNP_SUV_info_retriever(dicomdir)

% Loading dicom header into info
list=dir(dicomdir);
list=list(~startsWith({list.name}, '.')); % removing hidden files from the list
dicomimg=fullfile(dicomdir,list(10).name); % grabs tenth image (pseudorandom)
info = dicominfo(dicomimg,'UseDictionaryVR',true);

% Grab patient name and id
if isfield(info.PatientName,'FamilyName') && isfield(info.PatientName,'GivenName') 
    ptname=[info.PatientName.FamilyName '_' info.PatientName.GivenName];
else
    ptname='unknown, please refer to ptID';
end

if isfield(info,'PatientID')
    ptID=info.PatientID;
else
    ptID='unknown, please refer to ptname';
end

% Grab bodyweight
if isfield(info,'PatientWeight')
    bodyweight = 1000*info.PatientWeight; % in g
    if bodyweight <= 0 || bodyweight > 300000
        error('PatientWeight specified in dcm not plausible')
    end
else
    error('no PatientWeight specified in dcm')
end

% Grab injection dosis
injdosis=0;
while injdosis==0
    if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
        if info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose ~= 0
            injdosis = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % in Bq
        else
            info.RadiopharmaceuticalInformationSequence.Item_1 = rmfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose');
        end
    elseif isfield(info,'Private_0009_1038')
        if info.Private_0009_1038 ~= 0
            injdosis = 1000000*info.Private_0009_1038;
        else
            info=rmfield(info,'Private_0009_1038');
        end
    elseif isfield(info,'Private_0009_103a')
        if info.Private_0009_103a ~= 0
            injdosis = 1000000*info.Private_0009_103a;
        else
            info=rmfield(info,'Private_0009_103a');
        end
    elseif isfield(info,'RadionuclideTotalDose')
        if info.RadionuclideTotalDose ~= 0
            injdosis = info.RadionuclideTotalDose;
        else
            info=rmfield(info,'RadionuclideTotalDose');
        end
    else
        error('no RadionuclideTotalDose specified in dcm')
    end
end

% Get acquisition time
if isfield(info,'AcquisitionTime')
    acqtime= info.AcquisitionTime;
    if length(acqtime)>6
        acqtime=acqtime(1:6);
    end
else
    error('no AcquisitionTime specified in dcm')
end

% Get injection time
if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadiopharmaceuticalStartTime')
    injtime= info.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
    if length(injtime)==6
        injtime=[injtime '.000'];
    end
elseif isfield(info,'Private_0009_1039')
    injtime= info.Private_0009_1039;
    if length(injtime)==6
        injtime=[injtime '.000'];
    end
elseif isfield(info,'Private_0009_103b')
    injtime= info.Private_0009_103b;
    if length(injtime)==17
        injtime=[injtime(9:end) '0'];
    end
elseif isfield(info,'Private_0009_103d')
    injtime= info.Private_0009_103d;
    if length(injtime)==17
        injtime=[injtime(9:end) '.00'];
    end
else
    error('no RadiopharmaceuticalStartTime specified in dcm')
end

% Get halftime radionuclide
if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideHalfLife')
    halftime=info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife; % in seconds
elseif isfield(info,'RadionuclideHalfLife')
    halftime=info.RadionuclideHalfLife;
else
    error('no RadionuclideHalfLife specified in dcm')
end

% Get delta between acquistion time and injection time
dt=datetime(acqtime,'InputFormat','HHmmss')-datetime(injtime,'InputFormat','HHmmss.SSS');
dt=seconds(dt); % in seconds
if dt<1800
    error('difference between acquisition time and injection time is doubtfully small')
end      

% Calculate rescale factor
rescale_factor=(bodyweight/injdosis) * (1/2)^(-dt/halftime);

