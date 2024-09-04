  %This script runs through all neurons in all anesthetic conditions
%to determine their entrainment to the delta portion of the CSD
%5/23/23 JL

%% loading data
expTypes = {'Awake', 'LowIso', 'Ket'};

%load layer categorization
load('allLayerCats.mat')  

%initialize cell matrix for MI, p value, z-score for all
%2 timings (pre- and post- stim) by 3 conditions (awake, iso, ket)
allRealMIs = cell(2,3); %every cell's modulation index 
allPVals = cell(2,3); %every cell's p value of the MI compared to shuffled surrogates
allZScores = cell(2,3); %every cell's z score of the MI compared to shuffled surrogates

totalCells = zeros(1, 3); %counter for total cells 
%% determine entrainment 
%for each anestetic condition
for condition = 1%:3 
    %set data path
    dataPath = ['..\allData\All_Mouse_Data\',...
        expTypes{condition}, '\'];
    experimentFolders = dir([dataPath, '\CB*']);
    numExps = numel(experimentFolders);
    
    %prestim is row 1
    allRealMIs{1, condition} = cell(1,numExps);
    allPVals{1, condition} = cell(1,numExps);
    allZScores{1, condition} = cell(1,numExps);
    
    %poststim is row 2
    allRealMIs{2, condition} = cell(1,numExps);
    allPVals{2, condition} = cell(1,numExps);
    allZScores{2, condition} = cell(1,numExps);
    
    totalCells(condition) = 0;

    %for each experiment
    for mouse = 1%:numExps
        phyResName = dir([dataPath, experimentFolders(mouse).name, '\CB*_phyOutput.mat']).name;
        mouseName = phyResName(1:end-14);
        %get useful variables
        [trialSpiketrain, trueLoc, unitLayerCat, spiketrain, trialFireRate] = getUsefulVars_allData(mouseName, expTypes{condition});
        numUnits = numel(spiketrain);
        
        %get delta
        [csdDelta, allStartTimes] = getCSDDelta_allData(mouseName, expTypes{condition});
        
        %initialize statistic matrices for this experiment
        allRealMIs{1, condition}{mouse} = zeros(1,numUnits);
        allPVals{1, condition}{mouse} = zeros(1,numUnits);
        allZScores{1, condition}{mouse}= zeros(1,numUnits);
        allRealMIs{2, condition}{mouse} = zeros(1,numUnits);
        allPVals{2, condition}{mouse}= zeros(1,numUnits);
        allZScores{2, condition}{mouse}= zeros(1,numUnits);
        
        %for each unit
        for unit = 1:numel(spiketrain)
            %if it isn't a low firing rate cell
             if(trialFireRate(unit) > 0.5)
                %ensure channel isnt noisy
                currChan =  trueLoc(unit)-64;
                while(nnz(isnan(csdDelta(:, currChan, :))))
                    currChan = currChan - 1;
                end
                
                %get MI prestim
                [~,~,currMI(1)] = betterPhaseProb(csdDelta(:, currChan, :), trialSpiketrain(unit,:), false, [100 1000]);
                allRealMIs{1, condition}{mouse}(unit) = currMI(1);
                %get MI poststim
                [~,~,currMI(2)] = betterPhaseProb(csdDelta(:, currChan, :), trialSpiketrain(unit,:), false);
                allRealMIs{2, condition}{mouse}(unit) = currMI(2);
            
                %find poisson probability over a sliding window
                %window size
                windowSize = 200;
                %window shift
                windowShift = 1;
                trialLen = size(csdDelta, 1); %length of trial
                numTrials = size(trialSpiketrain,2);
                
                probs = zeros(1, (trialLen - windowSize)/windowShift);
                counter = 1;
                allSpikes = cell2mat(trialSpiketrain(unit, :));
                %for every window shift over entire trial
                for k = (1:windowShift:trialLen - windowSize) - 1000
                    % count total number of neurons in this window
                    numNeurons = numel(find(allSpikes > k & allSpikes < (k + windowSize)));
                    % divide by window, number of trials
                    probs(counter) = numNeurons/windowSize/numTrials;
                    counter = counter + 1;
                end
                
                %create surrogate experiments using the possion firing rate
                manyMIs = zeros(2, 100);
                %for each surrogate 
                for surrogate = 1:numTrials
                    %generate poisson trials
                    surrogateTrial = genPoissonTrials(probs, size(trialSpiketrain, 2), (1:numel(probs))-1000);
                    %get MIs
                    %prestim
                    [~,~,manyMIs(1, surrogate)] = betterPhaseProb(csdDelta(:, currChan, :), surrogateTrial, false, [100 1000]);
                    %poststim
                    [~,~,manyMIs(2, surrogate)] = betterPhaseProb(csdDelta(:, currChan, :), surrogateTrial, false);
                end
                
                %record pval, zscore
                for timing = 1:2 %prestim and poststim
                    allPVals{timing, condition}{mouse}(unit) = normcdf(currMI(timing),...
                        mean(manyMIs(timing,:)), std(manyMIs(timing,:)), 'upper');
                    allZScores{timing, condition}{mouse}(unit) = (allRealMIs{timing, condition}{mouse}(unit)...
                        - mean(manyMIs(timing,:)))/std(manyMIs(timing,:));
                end
                %counter for how many cells were run
                totalCells(condition) = totalCells(condition) + 1;
                
                disp(unit);
                
            else
                allPVals{1, condition}{mouse}(unit) = NaN;
                allZScores{1, condition}{mouse}(unit) = NaN;
                allRealMIs{1, condition}{mouse}(unit) = NaN;
                allPVals{2, condition}{mouse}(unit) = NaN;
                allZScores{2, condition}{mouse}(unit) = NaN;
                allRealMIs{2, condition}{mouse}(unit) = NaN;
            end %end if it's a high firing rate cell
        end %end for each unit
    end %end for each exp
end %end for each condition

%% get significance of each cell
significance = cell(2,3);
for condition = 1:3 %for each condition
    numExps = (numel(allPVals{1, condition}));
    significance{1, condition} = cell(1, numExps);
    significance{2, condition} = cell(1, numExps);
    for mouse = 1:numExps %for each mouse
        numUnits = numel(allPVals{1, condition}{mouse});
        significance{1, condition}{mouse} = zeros(1, numUnits);
        significance{2, condition}{mouse} = zeros(1, numUnits);
        for unit = 1:numUnits %for each unit
            for timeSeg = 1:2 %for prestim and poststim
                currPVal = allZScores{timeSeg, condition}{mouse}(unit);
                if isnan(currPVal)
                    significance{timeSeg, condition}{mouse}(unit) = nan;
                else
                    %compare the p value to the bonferroni corrected p
                    %value
                    cutoff = icdf('normal',[0.05/totalCells(condition), 1 - 0.05/totalCells(condition)],0,1);
                    significance{timeSeg, condition}{mouse}(unit) = currPVal > cutoff(2);
                end %end of checking for NaN
            end %end of for prestim and poststim
        end %end of for each unit
    end %end of for each mouse
end %end of for each condition
%

%% save data
save('..\dataFiles\poissonSurrogateAllConds.mat', 'allRealMIs', 'allPVals', 'allZScores', 'significance')
            