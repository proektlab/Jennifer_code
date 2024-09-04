%This creates a correlogram for all pairings of an entrained v1 and entrained
%ppa cell in all good awake, ket and low iso expereiments

expTypes = {'Awake', 'LowIso', 'Ket100'};

load('poissonSurrogateAllConds.mat', 'significance');
allSigs = significance(2, :);

load('allLayerCats.mat')
allLayerCats = {awakeLayerCats, lowIsoLayerCats, ketLayerCats};

%% get data, instantiate variables
histBins = -200:4:200;
numBins = numel(histBins) - 1;

stackedPrestims = cell(1,3);
stackedPoststims = cell(1,3);
sortedStackedPrestims = cell(1,3);
sortedStackedPoststims = cell(1,3);
prestimCounts = cell(1,3);
poststimCounts = cell(1,3);
mouseNumber = cell(1,3);
ppaNumber = cell(1,3);
v1Number = cell(1,3);

allLags = cell(2, 3);
%% for each condition
for condition = 1:3
    dataPath = ['Z:\code\Jennifer_code\SpikeSorting\usefulScripts\anesPaperFigures\allData\All_Mouse_Data\', expTypes{condition}];
    experimentFolders = dir([dataPath, '\CB*']);
    
    %all prestim lags
    stackedPrestims{condition} = zeros(0, numBins);
    %all poststim lags
    stackedPoststims{condition} = zeros(0, numBins);
    prestimCounts{condition} = [];
    poststimCounts{condition} = [];
    mouseNumber{condition} = [];
    ppaNumber{condition} = [];
    v1Number{condition} = [];
    
    allLags{1, condition} = {};
    %for each mouse
    for mouse = 1:numel(experimentFolders)
        phyResName = dir([dataPath, '\', experimentFolders(mouse).name, '\CB*_phyOutput.mat']).name;
        fullExpName = phyResName(1:end-14);
        [trialSpiketrain, trueLoc, unitLayerCat, spiketrain, trialFireRate] =...
            getUsefulVars_allData(fullExpName, expTypes{condition}, 'LED100');
        
        %obtain csd
        csdDelta = getCSDDelta_allData(fullExpName, expTypes{condition});
        
        %number of trials
        numTrials = size(trialSpiketrain, 2);
        
        %
        noNaNLayer = allLayerCats{condition}{mouse};
        noNaNLayer(isnan(noNaNLayer)) = 0;
        noNaNSig = allSigs{condition}{mouse};
        noNaNSig(isnan(noNaNSig)) = 0;
        sigPPA = find((noNaNLayer > 5) & noNaNSig);
        sigV1 = find((noNaNLayer <= 5) & noNaNSig);
        %for each significant ppa cell
        for ppaCell = 1:numel(sigPPA)
            %for each significant v1 cell
            for v1Cell = 1:numel(sigV1)
                %get all lags, which is v1 cell - ppa cell
                %so if v1 cell fires before ppa cell, the lag will be negative
                %and if v1 cell fires after ppa cell, the lag will be positive
                chosenV1 = sigV1(v1Cell);
                chosenPPA = sigPPA(ppaCell);
                
                currPrestimLags = findLagsWithinTrials(trialSpiketrain(chosenV1, :), trialSpiketrain(chosenPPA, :), [-500 0]);
                currPoststimLags = findLagsWithinTrials(trialSpiketrain(chosenV1, :), trialSpiketrain(chosenPPA, :), [0 500]);
                
                if(~isempty(currPrestimLags) && ~isempty(currPoststimLags))
                    stackedPrestims{condition} (end + 1, :) = histcounts(currPrestimLags, histBins, 'Normalization', 'probability');
                    stackedPoststims{condition} (end + 1,:) = histcounts(currPoststimLags, histBins, 'Normalization', 'probability');
                    
                    prestimCounts{condition}(end + 1) = numel(currPrestimLags);
                    poststimCounts{condition}(end + 1) = numel(currPoststimLags);
                    mouseNumber{condition}(end + 1) = str2num(experimentFolders(mouse).name(3:4));
                    ppaNumber{condition}(end + 1) = chosenPPA;
                    v1Number{condition}(end + 1) = chosenV1;
                    
                    allLags{1, condition}{end + 1} = currPrestimLags;
                    allLags{2, condition}{end + 1} = currPoststimLags;
                end
            end
        end
    end
end
%% look at spiketrain firing rates
allPPAPrestim = {[], [], []};
allV1Prestim =  {[], [], []};

allPPAPoststim =  {[], [], []};
allV1Poststim =  {[], [], []};

v1MouseInd=  {[], [], []};
ppaMouseInd=  {[], [], []};
v1NeuronInd =  {[], [], []};
ppaNeuronInd =  {[], [], []};

for condition = 1:3
    dataPath = ['Z:\code\Jennifer_code\SpikeSorting\usefulScripts\anesPaperFigures\allData\All_Mouse_Data\', expTypes{condition}];
    experimentFolders = dir([dataPath, '\CB*']);
    for mouse = 1:numel(experimentFolders)
        phyResName = dir([dataPath, '\', experimentFolders(mouse).name, '\CB*_phyOutput.mat']).name;
        fullExpName = phyResName(1:end-14);
        %load spiketrain
        [trialSpiketrain] = getUsefulVars_allData(fullExpName, expTypes{condition}, 'LED100');
        %for classification into region and significance
        noNaNLayer = allLayerCats{condition}{mouse};
        noNaNLayer(isnan(noNaNLayer)) = 0;
        noNaNSig = allSigs{condition}{mouse};
        noNaNSig(isnan(noNaNSig)) = 0;
        sigPPA = find((noNaNLayer > 5) & noNaNSig);
        sigV1 = find((noNaNLayer <= 5) & noNaNSig);
        
        numTrials = size(trialSpiketrain,2);
        %obtain prestim firing rate
        %ppa
        countPPA = cellfun(@(x) nnz(x<0&x>-500), trialSpiketrain(sigPPA,:));
        
        firingPPA = sum(countPPA, 2)/numTrials/.5;
        allPPAPrestim{condition} = [allPPAPrestim{condition}, firingPPA'];
        
        %v1
        countV1 = cellfun(@(x) nnz(x<0&x>-500), trialSpiketrain(sigV1,:));
        firingV1 = sum(countV1, 2)/numTrials/.5;
        allV1Prestim{condition} = [allV1Prestim{condition}, firingV1'];
        
        %obtain poststim firing rate
        %ppa
        countPPA = cellfun(@(x) nnz(x>0&x<500), trialSpiketrain(sigPPA,:));
        firingPPA = sum(countPPA, 2)/numTrials/.5;
        allPPAPoststim{condition} = [allPPAPoststim{condition}, firingPPA'];
        
        %v1
        countV1 = cellfun(@(x) nnz(x>0&x<500), trialSpiketrain(sigV1,:));
        firingV1 = sum(countV1, 2)/100/.5;
        allV1Poststim{condition} = [allV1Poststim{condition}, firingV1'];
        
        %note mouse and cell
        v1MouseInd{condition} = [v1MouseInd{condition}, mouse*ones(size(firingV1'))];
        ppaMouseInd{condition} = [ppaMouseInd{condition}, mouse*ones(size(firingPPA'))];
        v1NeuronInd{condition} = [v1NeuronInd{condition}, sigV1];
        ppaNeuronInd{condition} = [ppaNeuronInd{condition}, sigPPA];
        
    end
end

%histogram
figure
for condFig = 1:3
    subplot(3, 2, 2*condFig - 1);
    histogram(allV1Prestim{condFig}, 0:2:20);
    title([expTypes{condFig}, "V1"])
    
    subplot(3, 2, 2*condFig);
    histogram(allPPAPrestim{condFig}, 0:2:20);
    title([expTypes{condFig}, "PPA"])
end

%%
%re-sort
resortedLags = cell(1, 3);
removedCellsInd = cell(1,3);
newStackedPrestims= stackedPrestims;
newStackedPoststims = stackedPoststims;
isLow = false; %true if firing rate too low

%for each condition
for condInd = 1:3
    %remove bad looking prestims
    %for each cell combo
    for comboInd = numel(mouseNumber{condInd}):-1:1
        currV1 = find(v1MouseInd{condInd} == mouseNumber{condInd}(comboInd)...
            &  v1NeuronInd{condInd} ==v1Number{condInd}(comboInd));
        currPPA = find(ppaMouseInd{condInd} == mouseNumber{condInd}(comboInd)...
            &  ppaNeuronInd{condInd} ==ppaNumber{condInd}(comboInd));
        
        isLow = allV1Prestim{condInd}(currV1) < 0.5| allPPAPrestim{condInd}(currPPA) < 0.5;
        
        %check firing rates
        if isLow
            newStackedPrestims{condInd}(comboInd, :) = [];
            newStackedPoststims{condInd}(comboInd, :) = [];
            removedCellsInd{condInd} = [removedCellsInd{condInd}, comboInd];
            isLow = false;
        end
    end
    
    %get the max value and time of max value of each row
    [maxVal, maxInd] = max(-newStackedPrestims{condInd}+ newStackedPoststims{condInd}, [], 2);
    %sort max values
    [~, newOrder] = sort(maxVal);
    %re-sort
    resortedLags{condInd} = -newStackedPrestims{condInd}(newOrder, :)+ newStackedPoststims{condInd}(newOrder, :);
    resortedPrestim{condInd} = newStackedPrestims{condInd}(newOrder, :);
    resortedPoststim{condInd} = newStackedPoststims{condInd}(newOrder, :);
end
%% make figure post stim and pre stim

%for each condition
for condFig = 2:3
    %get the subtracted correl
    currFigImage = resortedPrestim{condFig};
    %subplot
    subplot(2, 2, condFig+1);
    %make the image
    imagesc(histBins(1:end-1), 1:size(currFigImage, 1),  currFigImage);
    
    %hold on
    hold on;
    %draw line at x=0
    xline(0, 'LineWidth', 2, 'LineStyle', '--');
    %colorbar
    %labels
    ylabel("V1 Unit ID")
    xlabel("Time Delay (ms)")
    
    colormap('winter');
    c = colorbar;
    c.Label.String = 'Change In Probability of V1 Spikes';
    caxis([0, 0.015])
    title([expTypes{condFig}, "Prestim"])
    
    currFigImage = resortedPoststim{condFig};
    %subplot
    subplot(2, 2, condFig-1);
    %make the image
    imagesc(histBins(1:end-1), 1:size(currFigImage, 1),  currFigImage);
    
    %hold on
    hold on;
    %draw line at x=0
    xline(0, 'LineWidth', 2, 'LineStyle', '--');
    %colorbar
    %labels
    ylabel("V1 Unit ID")
    xlabel("Time Delay (ms)")
    
    colormap('winter');
    c = colorbar;
    c.Label.String = 'Change In Probability of V1 Spikes';
    caxis([0, 0.015])
    title([expTypes{condFig}, "Poststim"])
end
