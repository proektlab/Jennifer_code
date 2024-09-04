% this adds corrected channel numbers, spiketrain by trials and layer info
% to postkilosorted .mat file
% a post post kilosort

%% loading things
clear

%load both chan maps
load('64ChNNGrid_2xH4CNT_remappedPlex_allConnected.mat', 'chanMap')
remappedCM = chanMap;
flipRemapped = [remappedCM(1:64), flip(remappedCM(65:96)), flip(remappedCM(97:end))];

load('64ChNNGrid_2xH4CNT_notRemappedPlex.mat', 'chanMap')
notRemappedCM = chanMap;
flipNotRemapped = [notRemappedCM(1:64), flip(notRemappedCM(65:96)), flip(notRemappedCM(97:end))];

%diff mapping: 21, 22, 24, 26, 27
%rightside up: CB22, 24, 26, 27, 37
isUpsideDown = ones(1, 11);
isDiff = zeros(1, 11);
isDiff(2:3) = [1, 1];
isUpsideDown([2, 3, 9]) = 0;

expType = 'VEPs';



%layer categories
layerCats = {1, [2, 3], 4, 5, 6, [9, 10, 11, 12], [13, 14]};

%experiment conditions
anesConds = {'Awake', 'LowIso', 'Ket100'};
%% make some files

%for each anesthetic condition
for cond = 1:3
    
    allPath = strcat('..\allData\All_Mouse_Data\', anesConds(cond), '\CB*');
    allPath = allPath{1};
    experimentFolders = dir(allPath);
    numExp = numel(experimentFolders);
    
    %for each file
    for i = 1:numExp
        dataPath = allPath(1:end-3);
        %load post kilo sort file
        cd(strcat(dataPath, experimentFolders(i).name)); %go to that data folder
        phyResName = dir([dataPath, experimentFolders(i).name, '\CB*_phyOutput.mat']).name;
        fullExpName = phyResName(1:end-14);
        load(phyResName) %load  datafiles
        
        %all firing rates
        allFireRates = cell(1, numel(experimentFolders));
        
        %load layer data
        [~, ~, ~, ~, ~, probeInfo] = getUsefulVars_allData(fullExpName, anesConds{cond});
        [~, allStartTimes] = getCSDDelta_allData(fullExpName, anesConds{cond});
        
        fragInfo = probeInfo.fragInfo;
        pr1_Layers = probeInfo.pr1_Layers;
        pr2_Layers = probeInfo.pr2_Layers;
        
        %make a spiketrain that's per trial
        numTrials = numel(allStartTimes);
        numUnits = numel(channel);
        trialSpiketrain = cell(numUnits, numTrials - 1);
        correctedSpiketrain = cell(1, numUnits);
        pseudoTrialSpiketrain = cell(numUnits, numTrials - 1);
        
        numV1 = 0;
        numPPA = 0;
        
        %locate fragments
        isFrag = ~isempty(fragInfo) && fragInfo.numFrags > 1; %determine if there's fragments
        if isFrag
            disp(['fixing fragmented ', num2str(i)])
            lengthFragMat = cell2mat(fragInfo.lengthFrags);
            fragTrialsStart = zeros(1,fragInfo.numFrags - 1);
            for fragment = 1:fragInfo.numFrags-1
                fragTrialsStart(fragment) = find(allStartTimes > fragInfo.startTimes(fragment+1), 1);
            end
        end
        
        
        %for every unit in spike train
        for j =  2:numTrials
            for k = 1:numUnits
                %find all spikes around the trial starting time, in ms
                trialSpikes = spiketrain{k}(spiketrain{k} > allStartTimes(j).*1000 - 1000 &...
                    spiketrain{k} < allStartTimes(j).*1000 + 2000);
                pseudoTrialSpikes = spiketrain{k}(spiketrain{k} > allStartTimes(j).*1000 - 1600 &...
                    spiketrain{k} < allStartTimes(j).*1000 - 100);
                %add these spiketimes to trialSpiketrain, adjusting for starting time
                if isFrag
                    adjustment = sum((j > fragTrialsStart) .* lengthFragMat);
                    trialSpiketrain{k, j-1} = (trialSpikes - allStartTimes(j).*1000 + adjustment*1000)';
                    pseudoTrialSpiketrain{k, j-1} = (pseudoTrialSpikes - allStartTimes(j).*1000 + adjustment*1000)';
                    
                    %adjust spiketrain for any fragments
                    correctedSpiketrain{k} = spiketrain{k};
                    for fragment = 1:fragInfo.numFrags-1
                        currFragTime = fragInfo.startTimes(fragment +1)*1000;
                         correctedSpiketrain{k}(correctedSpiketrain{k} > currFragTime)=...
                            correctedSpiketrain{k}(correctedSpiketrain{k} > currFragTime)...
                            - fragInfo.lengthFrags{fragment+1}*1000;
                    end
                    
                else
                    trialSpiketrain{k, j-1} = (trialSpikes - allStartTimes(j).*1000)';
                    pseudoTrialSpiketrain{k, j-1} = (pseudoTrialSpikes - allStartTimes(j).*1000)';
                    correctedSpiketrain{k} = spiketrain{k};
                end
            end
        end
        
        clear fragInfo
        
        %average firing rate in trials?
        trialFireRate = sum(cellfun(@numel, trialSpiketrain), 2) ./ (3*(numTrials-1));
        allFireRates{i} = trialFireRate;
        
        %determine correct chanmap
        if(isDiff(i)) %if is different mapping
            if(isUpsideDown(i)) %and upside down
                currCM = flipNotRemapped;
            else %not upside down
                currCM = notRemappedCM;
            end
        else %same mapping
            if(isUpsideDown(i)) %and upsidedown
                currCM = flipRemapped;
            else %not upside down
                currCM = remappedCM;
            end
        end
        
        
        trueLoc = zeros(1, numUnits);
        %something like this:
        for j = 1:numUnits
            trueLoc(j) = find(currCM == channel(j));
        end
        
        
        %alter layer assign to affect all of the probe
        %if this mouse has layer data
        pr1Assigned = true;
        pr2Assigned= true;
        if ~exist('pr1_Layers', 'var')
            pr1_Layers.LayerAssign = zeros(32, 1);
            pr1Assigned = false;
        end
        if ~exist('pr2_Layers', 'var')
            pr2_Layers.LayerAssign = zeros(32, 1);
            pr2Assigned = false;
        end
        %put layer assigns together
        bigLayerAssign = [pr1_Layers.LayerAssign', pr2_Layers.LayerAssign' + 8];
        unitLayers = zeros(1, numUnits);
        %for each unit
        for j = 1:numUnits
            %determine layer with layerassign
            unitLayers(j) = bigLayerAssign(trueLoc(j)-64);
        end
        
        %determine layer category of each unit?
        %for each unit
        unitLayerCat = zeros(1, numUnits);
        for j = 1:numUnits
            %get layer
            currLayer = find(cell2mat(cellfun(@(x) any(x == unitLayers(j)), layerCats, 'UniformOutput', false)));
            %empty means its out of cortex, also remove cells with low firerate
            if (isempty(currLayer) || trialFireRate(j) < 0.5)
                unitLayerCat(j) = NaN;
            else
                unitLayerCat(j) = currLayer;
                if currLayer < 6
                    numV1 = numV1 + 1;
                else
                    numPPA = numPPA + 1;
                end
            end
        end
        
        isGood = ones(1, numUnits);
        %determine if good (in cortex and high firing rate)
        for j = 1:numUnits
            if(isnan(unitLayerCat(j)) & trialFireRate(j) < 0.5)
                isGood(j) = 0;
            end
        end
        
        %make save name
        saveName = [fullExpName, '_usefulVars'];
        %(make sure to sancheck first) save everything back into file
        save(saveName, 'trialSpiketrain', 'pseudoTrialSpiketrain', 'trialFireRate', 'unitLayers', 'trueLoc', 'unitLayerCat', 'isGood', 'numV1', 'numPPA', 'correctedSpiketrain', '-append');
    end
    
end
