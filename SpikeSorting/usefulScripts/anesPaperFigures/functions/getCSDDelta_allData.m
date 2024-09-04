%function [csdDelta, allStartTimes] = getCSDDelta(mouseName, expType, stimType)
%loads CSD delta of a mouse in iso_awake_VEPs
%inputs: mouseName: name of mouse (eg 'CB34_Awake_LED100')
%        expType: category of experiment ('Awake', 'LowIso' or 'Ket')
%        stimType:  type of stimulus ('LED100', 'FSCF' or 'Odd')
%outputs: csdDelta: delta filtered CSD. timepoints x channel x trials
%         allStartTimes: all trial start times

%jl 6/29/23
function [csdDelta, allStartTimes, csdGamma] = getCSDDelta_allData(mouseName, expType)
csdDelta = [];
allStartTimes = [];
fileLoc = ['..\allData\All_Mouse_Data\',...
    expType, '\', mouseName(1:4), '\'];
%go to correct layer data location
if exist(fileLoc, 'file')
    load([fileLoc, mouseName, '_CSD.mat'], 'allStartTimes', 'csdDelta', 'csdGamma') %load trialtimes
else
    disp("File not found!")
end
end