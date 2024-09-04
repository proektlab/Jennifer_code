%function [trialSpiketrain, trueLoc, unitLayerCat, spiketrain, trialFireRate] =...
%    getUsefulVars(mouseName, expType, stimType)
%function that loads spike data about an experiment from its files
%inputs:    mouseName: name of mouse (eg 'CB34_Awake_LED100')
%           expType: category of experiment ('Awake', 'LowIso' or 'Ket')
%                   an expType of 'oddballMerge' accesses the merged
%                   oddball experiments which have awake and iso
%           stimType: type of stim, 'LED' or 'FSCF'
%output:    trialSpiketrain: units x trials cell containing spike times in ms in
%           each trial, relative to stimulus.
%           trueLoc: vector containing true channel of each cell
%           unitLayerCat: vector containing layer category of each unit. 1
%           is v1 layer 1, 2 is v1 layer 2/3, 3 is v1 layer4, 4 is v1 layer
%           5, 5 is v1 layer 6, 6 is upper ppa, 7 is lower ppa
%           spiketrain: times of spiking for every unit over the entire
%           experiment
%           trialFireRate: firing rate of each trial
%           probeData: variable containing data on recording, including the
%           layer categories of the probe during that experiment
%jl 8/5/21
function [trialSpiketrain, trueLoc, unitLayerCat, spiketrain, trialFireRate, probeData, spikeData] =...
    getUsefulVars_allData(mouseName, expType, stimType)
if(nargin < 3)
    stimType = 'LED100';
end
trialSpiketrain = [];
unitLayerCat = [];
spiketrain = [];
trueLoc = [];
trialFireRate = [];
probeData = [];

%allow for acceptance of both 'Ket100' and 'Ket' to refer to ketamine experiments 
if strcmp(expType, "Ket")
    expType = "Ket100";
end

fileLoc = strcat("..\allData\All_Mouse_Data\",...
    expType, "\", mouseName(1:4), "\");



%check for file
if (exist(fileLoc, 'dir'))
    load(strcat(fileLoc, mouseName(1:4), "_", expType, "_", stimType, "_phyOutput.mat"), 'spiketrain');
    load(strcat(fileLoc, mouseName(1:4), "_", expType, "_", stimType, "_usefulVars.mat"),...
        'unitLayerCat', 'trialSpiketrain', 'trueLoc', 'trialFireRate', 'isGood', 'meanWaveform', 'correctedSpiketrain')
    load(strcat(fileLoc, mouseName(1:4), "_", expType, "_", stimType,...
        "_ProbeData.mat"), 'pr1_Layers', 'pr2_Layers', 'fragInfo');
    probeData.pr1_Layers = pr1_Layers;
    probeData.pr2_Layers = pr2_Layers;
    spikeData.meanWaveform = meanWaveform;
    spikeData.isGood = isGood;
    spikeData.correctedSpiketrain = correctedSpiketrain;
    if(exist('fragInfo', 'var'))
        probeData.fragInfo = fragInfo;
    else
        probeData.fragInfo = [];
    end
else
    disp("File not found!")
end
end