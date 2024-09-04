%This script finds if a unit has a significant increase or decrease in
%firing rate after the stimulus
%Then, it generates a figure showing the proportion of units in all
%conditions that change firing rate and determines if the proportion of
%units that change firing rate is significantly different across anesthetic
%conditions.
%JL 5/23/2023

%% Determine if each unit changes in firing rate
% do for each condition
expTypes = {'Awake', 'LowIso', 'Ket100'};
load poissonSurrogateAllConds.mat
totCells = [132, 201, 50; 132, 231, 95]; %total cells per condition

%total cells for bonferroni correction
windowSize = 500;
sigFiringRate = nan(2, 3, 231);
firingZScores = nan(2, 3, 231);
mouseCounter = nan(2, 3, 231);
mouseNum = 0;
cellNum = nan(2, 3, 231);

cellCounter = zeros(2, 3); %for counting cells

%for each anesthetic condition
for condInd = 1:3
    %navigate to data
    dataPath = ['Z:\code\Jennifer_code\SpikeSorting\usefulScripts\anesPaperFigures\allData\All_Mouse_Data\', expTypes{condInd}];
    experimentFolders = dir([dataPath, '\CB*']);
    numExp = numel(experimentFolders);
    
    %for each experiment
    for mouse = 1:numExp
        mouseNum = mouseNum + 1;
        %load things
        phyResName = dir([dataPath, '\', experimentFolders(mouse).name, '\CB*_phyOutput.mat']).name;
        mouseName = phyResName(1:end-14);
        %get useful variables
        
        
        [trialSpiketrain, trueLoc, unitLayerCat, spiketrain, trialFireRate, probeData] = getUsefulVars_allData(mouseName, expTypes{condInd});
        numUnits = numel(spiketrain);
        %possible locations of all start times
        [~, allStartTimes] = getCSDDelta_allData(mouseName, expTypes{condInd});
        
        numTrials = numel(allStartTimes) - 1;
        %for each unit
        for unit = 1:numUnits
            %if it has a reasonable firing rate
            if(~isnan(unitLayerCat(unit)))
                realFRs = zeros(1,numTrials);
                surrogateMeans = zeros(1,100);
                currSpiketrain = spiketrain{unit};
                %determine probe
                if unitLayerCat(unit) < 6
                    probeInd = 1; %v1
                else
                    if unitLayerCat(unit) >= 6
                        probeInd = 2; %ppa
                    else
                        continue
                    end
                end
                %add to counter
                cellCounter(probeInd, condInd) = cellCounter(probeInd, condInd) + 1;
                %adjust for fragments
                isFrag = ~isempty(probeData.fragInfo); %determine if there's fragments
                if isFrag
                    if (probeData.fragInfo.numFrags) > 1
                        disp('fixing frag')
                        for i = 1:probeData.fragInfo.numFrags-1 %for each frag
                            currSpiketrain(currSpiketrain > (probeData.fragInfo.startTimes(i+1)*1000)) =...
                                currSpiketrain(currSpiketrain > (probeData.fragInfo.startTimes(i+1)*1000)) + probeData.fragInfo.lengthFrags{i+1};
                        end
                    end
                end
                probeData.fragInfo = [];
                
                
                %for each real trial
                for realTrial = 1:numTrials
                    currWindow = [allStartTimes(realTrial + 1)*1000,...
                        allStartTimes(realTrial + 1)*1000+windowSize];
                    %get the firing rate in this window
                    realFRs(realTrial) = nnz(currSpiketrain > currWindow(1) &...
                        currSpiketrain < currWindow(2))/windowSize;
                end
                %bootstrap the real trials
                bootstrapMeans = bootstrp(100, @mean, realFRs);
                %for 100 times
                for i = 1:100
                    surrogateFRs = zeros(1,100);
                    %for all fake trials
                    for j = 1:numTrials
                        randTime = rand * max(currSpiketrain)-windowSize;
                        currWindow = [randTime, randTime + windowSize];
                        %get the firing rate in this window
                        surrogateFRs(j) = nnz(currSpiketrain > currWindow(1) &...
                            currSpiketrain < currWindow(2))/windowSize;
                    end
                    surrogateMeans(i) = mean(surrogateFRs);
                end
                %t test
                [~, p] = ttest2(bootstrapMeans, surrogateMeans,  'Tail','right');
                %bonferroni correction
                sigFiringRate(probeInd, condInd, cellCounter(probeInd, condInd)) = p < 0.05/sum(totCells(:, condInd));
                [~, p2] = ttest2(bootstrapMeans, surrogateMeans, 'Tail' ,'left');
                if (p2 < 0.05/totCells(condInd)) %check for sig decreases
                    sigFiringRate(probeInd, condInd, cellCounter(probeInd, condInd))  = -1;
                end
                
                %record z value
                firingZScores(probeInd, condInd, cellCounter(probeInd, condInd)) =  -icdf('normal', p, 0, 1);
                mouseCounter(probeInd, condInd, cellCounter(probeInd, condInd)) = mouseNum;
                
            %else
                %sigFiringRate(probeInd, condInd, cellCounter(probeInd, condInd))  = NaN;
                %firingZScores(probeInd, condInd, cellCounter(probeInd, condInd)) = NaN;
            end
            cellNum(probeInd, condInd, cellCounter(probeInd, condInd)) = unit;
        end
    end
end

firingZScores(isinf(firingZScores)) = NaN; %replace infs with nans

%% save some variables in table format
FiringZScore = cell(32, 1);
SigFiring = cell(32, 1);
%for each mouse
for mouse = 1:32
   SigFiring{mouse} = sigFiringRate(mouseCounter == mouse);
   FiringZScore{mouse} = firingZScores(mouseCounter == mouse);
   [~, ind] = sort(cellNum(mouseCounter == mouse));
   
   FiringZScore{mouse} = FiringZScore{mouse}(ind);
   SigFiring{mouse} = SigFiring{mouse}(ind);
   SigFiring{mouse}(isnan(FiringZScore{mouse})) = NaN;
   

end
firing_stats_table = table(FiringZScore, SigFiring);
save("firing_stats_table.mat", "firing_stats_table");

%% bar plot
firingCount = {zeros(3, 3), zeros(3,3)}; %raw count
firingCountNorm = {zeros(3, 3), zeros(3,3)}; %fraction
firingCountUpper = {zeros(3, 3), zeros(3,3)}; %error bar from bootstrapping
firingCountLower = {zeros(3, 3), zeros(3,3)}; %error bar from bootstrapping

titles = {"V1", "PPA"}
firingTypeVal = [1 0 -1]; %increase, no change, decrease

fracFunc = @(x, numCells) nnz(x)/numCells;
for i = 1:2 %for each probe
    for j = 1:3 %for each condition
        %for each firing count
        for firingType = 1:3
            %count number of increasing firing rate
            firingCount{i}(j, firingType) = nnz(sigFiringRate(i, j, :) == firingTypeVal(firingType));
            %bootstrap it
            bootRes = bootci(100, @(x) fracFunc(x, totCells(i, j)) , sigFiringRate(i, j, :) == firingTypeVal(firingType));
            firingCountLower{i}(j, firingType) = bootRes(1);
            firingCountUpper{i}(j, firingType) = bootRes(2);
        end
        firingCountNorm{i}(j, :) = firingCount{i}(j, :)./totCells(i, j);
        
    end
end

allBar = cell(1, 2);
ff = figure('Color', 'w', 'Renderer', 'Painters');
for i = 1:2
    subplot(1, 2, i) 
    % Calculate the number of groups and number of bars in each group
    allBar{i} = bar(firingCountNorm{i}');
    hold on
    [ngroups,nbars] = size(firingCountNorm{i});
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for j = 1:nbars
        x(j,:) = allBar{i}(j).XEndPoints;
    end
    
    %bar graph number of firing rate changes
    firingError{i} = errorbar(x', firingCountNorm{i}', firingCountNorm{i}'...
        - firingCountLower{i}', firingCountUpper{i}' - firingCountNorm{i}', '.k');
    allBar{i}(1).FaceColor = [0.6350 0.0780 0.1840]; %red
    allBar{i}(2).FaceColor = [0.4940 0.1840 0.5560]; %purple
    allBar{i}(3).FaceColor = [0.4660 0.6740 0.1880]; %dark green
    %labels
    xticklabels(["Increase", "No Change", "Decrease"])
    title(titles(i));
    xlabel("Condition")
    ylabel("Fraction of Units")
    ylim([0, 1])
end

%legend
legend({"Awake", "Iso", "Ketamine"});
%title


%% statistical analysis of proportion of each condition that increases firing rate
numInc = zeros(2, 3);
chi2stat_inc = zeros(1, 2);
p_inc = zeros(1, 2);
%for each probe
for probeInd = 1:2
    %for each condy
    for condInd = 1:3
        %get the number of 1s and -1s
        numInc(probeInd, condInd) = nnz(sigFiringRate(probeInd, condInd, :) == 1);
    end
    [chi2stat_inc(probeInd),p_inc(probeInd)] = chiSquareTest([numInc(probeInd, :); cellCounter(probeInd, :)]);
    tmcomptest([numInc(probeInd, :); cellCounter(probeInd, :)]', 0.05)
end


%% statistical analysis of proportion of each condition that decreases firing rate
numDec = zeros(2, 3);
chi2stat_dec = zeros(1, 2);
p_dec = zeros(1, 2);
%for each probe
for probeInd = 1:2
    %for each condy
    for condInd = 1:3
        %get the number of 1s and -1s
        numDec(probeInd, condInd) = nnz(sigFiringRate(probeInd, condInd, :) == -1);
    end
    [chi2stat_dec(probeInd),p_dec(probeInd)] = chiSquareTest([numDec(probeInd, :); cellCounter(probeInd, :)]);
    tmcomptest([numDec(probeInd, :); cellCounter(probeInd, :)]', 0.05)
end

%% statistical analysis of proportion of each condition that changes firing rate
numChange = zeros(2, 3);
chi2stat_change = zeros(1, 2);
p_change = zeros(1, 2);
%for each probe
for probeInd = 1:2
    %for each condy
    for condInd = 1:3
        %get the number of 1s and -1s
        numChange(probeInd, condInd) = nnz(sigFiringRate(probeInd, condInd, :) == 0);
    end
    [chi2stat_change(probeInd),p_change(probeInd)] = chiSquareTest([numChange(probeInd, :); cellCounter(probeInd, :)]);
    tmcomptest([numDec(probeInd, :); cellCounter(probeInd, :)]', 0.05)
end

%% 
[chi2stat_total,p_total] = chiSquareTest([sum(numChange); sum(cellCounter)])