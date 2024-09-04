%This script counts the number of entrained units in each condition
%Then, it generates a figure showing the proportion of units in all
%conditions that are entrained and determines if the proportion of
%units that are entrained is significantly different across anesthetic
%conditions.
%JL 5/23/2023


%% load data and intialize matrices 
load('poissonSurrogateAllConds.mat', 'significance', 'allPVals', 'allZScores')

load('allLayerCats.mat')
allLayerCats = {awakeLayerCats, lowIsoLayerCats, ketLayerCats};

%matrix for barplot: condition by entrained
allBarData = zeros(3, 2);
v1Data = zeros(3,2);
ppaData = zeros(3,2);

%matrices for errorbars
ppaLow = zeros(3, 2);
ppaHigh = zeros(3, 2);
v1Low = zeros(3, 2);
v1High = zeros(3, 2);

totalSpikes = zeros(1,3);
totalV1 = zeros(1,3);
totalPPA = zeros(1,3);

thingCounter = 0; %for every condition and times

bootDists = zeros(12, 100);
allSigs = cell(1, 12);
allZs = cell(1,12);

sigsCombinedCount = zeros(1, 3); %for number of cells that are significant pre OR post
sigsCombined = cell(1,3);

%% count number of entrained cells in each probe
%for each condition
for condition = 1:3
    currBarData = [0 0];
    currV1Data = [0 0];
    currPPAData = [0 0];
    
    %for botstrapping
    currV1Sigs = cell(1, 2);
    currV1Zs = cell(1,2);
    currPPASigs = cell(1, 2);
    currPPAZs = cell(1,2);
    allLayNum{condition} = zeros(1,7);
    
    currCombinedSigsCount = 0;
    currCombinedSigs = [];
    %for each mouse
    for mouse = 1:numel(significance{2,condition})
        
        noNaNSig = significance{2, condition}{mouse}; %remove nans from significance
        noNaNSig(isnan(noNaNSig)) = 0;
        
        goodSigs = significance{2, condition}{mouse}; %for counting total 
        goodSigs(isnan(goodSigs)) = [];
        
        goodSigsPre = significance{2, condition}{mouse};
        goodSigsPre(isnan(goodSigsPre)) = [];
        
        
        %get v1 and ppa units
        v1Units = find(allLayerCats{condition}{mouse} < 6 & ~isnan(allPVals{2, condition}{mouse}));
        ppaUnits = find(allLayerCats{condition}{mouse} >= 6 & ~isnan(allPVals{2, condition}{mouse}));
        
        for layer = 1:7
            currLayNum(layer) = nnz(allLayerCats{condition}{mouse} == layer & noNaNSig);
        end
        allLayNum{condition} = allLayNum{condition} + currLayNum;
        %count up units
        totalSpikes(condition) = totalSpikes(condition) + numel(~isnan(allLayerCats{condition}{mouse} ));
        totalV1(condition) = totalV1(condition)+ numel(v1Units);
        totalPPA(condition) = totalPPA(condition) + numel(ppaUnits);
        for timeLoc = 1:2 %for prestim and poststim
            currBarData(timeLoc) = currBarData(timeLoc) + nnz(significance{timeLoc,condition}{mouse} == 1);
            %do that for v1
            currV1Data(timeLoc) = currV1Data(timeLoc) + nnz(significance{timeLoc,condition}{mouse}(v1Units)==1);
            %obtain significances for bootstrappin
            currV1Sigs{timeLoc} = [currV1Sigs{timeLoc}, significance{timeLoc,condition}{mouse}(v1Units)==1];
            currV1Zs{timeLoc} = [currV1Zs{timeLoc}, allZScores{timeLoc,condition}{mouse}(v1Units)];
            %do that for ppa
            currPPAData(timeLoc) = currPPAData(timeLoc) + nnz(significance{timeLoc,condition}{mouse}(ppaUnits)==1);
            currPPASigs{timeLoc} = [currPPASigs{timeLoc}, significance{timeLoc,condition}{mouse}(ppaUnits)==1];
            currPPAZs{timeLoc} = [currPPAZs{timeLoc}, allZScores{timeLoc,condition}{mouse}(ppaUnits)];
        end
        currCombinedSigsCount = currCombinedSigsCount + nnz(goodSigs | goodSigsPre);
        currCombinedSigs = [currCombinedSigs, goodSigs | goodSigsPre];
        
    end
    allBarData(condition, :) = currBarData/totalSpikes(condition);
    
    
    allV1Data(condition, :) = currV1Data/totalV1(condition);
    allPPAData(condition, :) = currPPAData/totalPPA(condition);
    
    
    %place all conditions into a large cell matrix for stats
    for timeLoc = 1:2
        thingCounter = thingCounter + 1;
        
        [bootRes, bootDists(thingCounter, :)] = bootci(100, @nanmean, currV1Sigs{timeLoc} );
        v1Low(condition, timeLoc) = bootRes(1);
        v1High(condition, timeLoc) = bootRes(2);
        allSigs{thingCounter} = currV1Sigs{timeLoc};
        allZs{thingCounter} = currV1Zs{timeLoc};
        
        thingCounter = thingCounter + 1;
        
        [bootRes,  bootDists(thingCounter, :)] = bootci(100, @nanmean,  currPPASigs{timeLoc});
        ppaLow(condition, timeLoc) = bootRes(1);
        ppaHigh(condition, timeLoc) = bootRes(2);
        allSigs{thingCounter} = currPPASigs{timeLoc};
        allZs{thingCounter} = currPPAZs{timeLoc};
    end
    sigsCombinedCount(condition) = currCombinedSigsCount;
    sigsCombined{condition} = currCombinedSigs;
end

distNames = ["awake pre v1", "awake pre ppa", "awake post v1", "awake post ppa",...
    "iso pre v1", "iso pre ppa", "iso post v1", "iso post ppa",...
    "ket pre v1", "ket pre ppa", "ket post v1", "ket post ppa"];
%% create bar plot with fraction of units over each condition
ff = figure('Color', 'w', 'Renderer', 'Painters');
subplot(1,2, 1)
allBar = bar(allV1Data);
allBar(2).FaceColor = [.2 .6 .5]; %green
allBar(1).FaceColor = [0.5 0.5 0.5]; %gray
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(allV1Data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = allBar(i).XEndPoints;
end

errorbar(x', allV1Data, allV1Data - v1Low, v1High - allV1Data, '.k')
%legend(["Prestim","Poststim"]);
xlabel("Condition")
xticklabels(["Awake", "Low Iso", "Ket"])
ylabel("Fraction of Units")
title('V1')

%show ppa units
subplot(1,2,2)
allBar = bar(allPPAData);
allBar(2).FaceColor = [.2 .6 .5]; %green
allBar(1).FaceColor = [0.5 0.5 0.5]; %gray


x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = allBar(i).XEndPoints;
end
hold on
errorbar(x', allPPAData, allPPAData - ppaLow, ppaHigh - allPPAData, '.k')

legend(["Prestim","Poststim"]);
xlabel("Condition")
xticklabels(["Awake", "Low Iso", "Ket"])
ylabel("Fraction of Units")
ylim([0 0.6])
title('PPA')

sgtitle('Entrainment of Units by Probe')


%% perform stats to determine if fraction of entrained units is significant 

%adjust allSigs so it fits the kruskalwallis function
kruskalSigs = cell2mat(cellfun(@(x) horzcat(x, NaN(1, 231 - numel(x)))', allSigs, 'UniformOutput', false));

tukeyVals = zeros(12, 2);
for i = 1:12
    tukeyVals(i, 1) = nnz(kruskalSigs(:, i) == 1);
    tukeyVals(i, 2) = nnz(~isnan(kruskalSigs(:, i)));
end

V1Cols = [1,3,5,7,9,11];
PPACols = [2,4,6,8,10,12];
tmcomptest(tukeyVals(V1Cols, :))

tmcomptest(tukeyVals(PPACols, :))

Vals4ChiSqr(1,:) = tukeyVals(:,1)';
Vals4ChiSqr(2,:) = (tukeyVals(:,2)- tukeyVals(:,1))';

%chi square test

[chi2stat_PPA,p_PPA] = chiSquareTest(Vals4ChiSqr(:,PPACols))

[chi2stat_V1,p_V1] = chiSquareTest(Vals4ChiSqr(:,V1Cols))

%% perform chi square  on total number of cells
[chi2stat_All,p_All] = chiSquareTest([sigsCombinedCount; totalSpikes])
%% perform tukey 
tmcomptest([sigsCombinedCount', totalSpikes']);