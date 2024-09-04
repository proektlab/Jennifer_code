%do some across probe stats? idk
clear
dataPath = 'D:\adeeti\KiloSort_results\LEDAwakeOnly\'; %CHANGE THIS
experimentFolders = dir([dataPath, 'Output*']);
layerCats = [1 2 4 5 6];
numLayers = numel(layerCats);
numExps = numel(experimentFolders);
thresh = 4;

load("poissonSurrogateDataAwakeLEDDiffProbe.mat", 'zScores');
%the matrix: 5 x 5 x 2 cell,
% rows are where neuron is from, columns is where the csd is from, z is v1
% or ppa
zScoresByLayer = cell(5, 5, 2);
numberCells = zeros(5, 2, numExps);

%for each mouse
for mouse = 1:numExps
    %load layer categories
    load([dataPath, experimentFolders(mouse).name, '\Data\usefulVars.mat'], 'unitLayers', 'trialFireRate')
    %for neuron
    for unit = 1:numel(unitLayers)
        %get layer
        %determine if v1 or ppa
        if any(unitLayers(unit) == [0 7 8 15]) %if not in cortex, break
            continue
        else %in cortex
            if unitLayers(unit) < 8
                probeCat = 1; %v1
                currUnitLay = unitLayers(unit);
                
            else
                probeCat = 2; %ppa
                currUnitLay = unitLayers(unit) - 8;
            end
        end
        
        currUnitLay = find(layerCats == currUnitLay); %because 2/3 combined :/
        
        %for each layer of csd (1, 2, 4, 5, 6) it was matched with
        for csdLayer = 1:numLayers
            %add the z score to the appropriate part of the matrix
            zScoresByLayer{csdLayer,currUnitLay, probeCat} =...
                [zScoresByLayer{csdLayer,currUnitLay, probeCat}, zScores{mouse}(layerCats(csdLayer), unit)];
            numberCells(currUnitLay, probeCat, mouse) = numberCells(csdLayer, probeCat, mouse)+1;
        end
        
    end
end

%% figure?
meanZ = cellfun(@(x) nnz(x > thresh)/nnz(~isnan(x)), zScoresByLayer);
%maybe a color plot for the mean of the cells? idk
figure('color', 'white', 'renderer', 'painters')
titles = {'V1 Cells w/ PPA CSD', 'PPA Cells w/ V1 CSD'};
for i = 1:2
    subplot(1, 2, i)
    currData = nanmean(meanZ(:,:,i));
    %currData(isnan(currData)) = 0;
    bar(1:5, currData);
    hold on
    %textlabels
    %text(b.XEndPoints,b.YEndPoints,string(currData),'HorizontalAlignment','center',...
    %'VerticalAlignment','bottom')
    title(titles{i})
    xlabel('Unit Layer')
    ylabel(['Fraction of Z Score Above ', num2str(thresh)]);
    xticks(1:5);
    xticklabels({'1', '2/3', '4', '5', '6'});
    %yticks(1:5);
    %yticklabels({'1', '2/3', '4', '5', '6'});
    %c = colorbar;
    %c.Label.String = ['Fraction of Z Score Above ', num2str(thresh)];
end

%legend(experimentFolders.name);
sgtitle('Cross Probe Entrainment');
%% figure but with e r r o r b a r

thresholdCross = cellfun(@(x) (rmmissing(x) > thresh), zScoresByLayer, 'UniformOutput', false);
low95 = zeros(2, 5);
high95 = zeros(2, 5);
%maybe a color plot for the mean of the cells? idk
figure('color', 'white', 'renderer', 'painters')
titles = {'V1 Cells w/ PPA CSD', 'PPA Cells w/ V1 CSD'};
allBootDists = zeros(2, 5, 100);

axesVals = {[0 0.65], [0, 0.3]};
for i = 1:2
    %for each layer (of the unit origin)
    for layInd = 1:5
        %get all the data from the layers across
        layData = cell2mat(thresholdCross(:, layInd, i)');
        %bootstrap
        if (~isempty(layData))
            [bootRes, allBootDists(i, layInd, :)] = bootci(100, @nanmean, layData);
            low95(i, layInd) = bootRes(1);
            high95(i, layInd) = bootRes(2);
        end
    end
    subplot(1, 2, i)
    currData = nanmean(meanZ(:,:,i));
    %currData(isnan(currData)) = 0;
    bar(1:5, currData);
    hold on
    errorbar(1:5, currData, currData - low95(i, :), currData - high95(i, :), '.k');
    hold on
    %textlabels
    %text(b.XEndPoints,b.YEndPoints,string(currData),'HorizontalAlignment','center',...
    %'VerticalAlignment','bottom')
    title(titles{i})
    xlabel('Unit Layer')
    ylabel(['Fraction of Z Score Above ', num2str(thresh)]);
    xticks(1:5);
    xticklabels({'1', '2/3', '4', '5', '6'});
    ylim(axesVals{i});
    %yticks(1:5);
    %yticklabels({'1', '2/3', '4', '5', '6'});
    %c = colorbar;
    %c.Label.String = ['Fraction of Z Score Above ', num2str(thresh)];
end

%legend(experimentFolders.name);
sgtitle('Cross Probe Entrainment');

%% to excel
figureSup16V1 = [nanmean(meanZ(:,:,1)); low95(1, :); high95(1, :)];
figureSup16PPA = [nanmean(meanZ(:,:,2)); low95(2, :); high95(2, :)];

figureSup16 = array2table([figureSup16V1,figureSup16PPA]);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figureSup16CrossProbe.xlsx';
writetable(figureSup16,filename,'Sheet',1,'Range','A1')
%%

for probe = 1:2
    tmInput = zeros(4, 2);
    for csdLayer = 2:numLayers
        currTMRow = csdLayer - 1;
        currZs = cell2mat(zScoresByLayer(:, csdLayer, probe));
        tmInput(currTMRow, 1) = nnz(currZs > thresh);
        tmInput(currTMRow, 2) = nnz(~isnan(currZs));
    end
    disp(titles(probe));
    tmcomptest(tmInput, .05)
end

%% t test the boot dists?
tVals = nan(2, 5, 5);
p_val = nan(2, 5, 5);
stds = nan(2, 5);

for probe = 1:2
    for cond1 = 1:5
        stds(probe, cond1) = nanstd(allBootDists(probe, cond1, :));
    end
    for cond1 = 1:5
        %for every other contrast
        for cond2 = cond1+1:5
            mean1 = nanmean(meanZ(:,cond1,probe));
            mean2 = nanmean(meanZ(:,cond2,probe));
            dist2 = allBootDists(probe, cond2, :); 


            %tVals(cond1, cond2) = (nanmean(dist1) - nanmean(dist2))/sqrt((nanstd(dist1)^2)/100 + (nanstd(dist2)^2)/100);
            %[~, ~, ~, tValsMatlab(cond1, cond2)] = ttest2(dist1', dist2');
            %p_val(cond1,cond2) = normcdf(mean1, nanmean(dist2), nanstd(dist2));
            tVals(probe, cond1, cond2) = (mean2 - mean1)/stds(probe, cond2);
            %disp((mean2 - mean1)/stds(probe, cond2));
        end
    end
end

disp(tVals);
pVals = (1 - tcdf(abs(tVals), 99))*10*2;

for pr = 1:2
    for i = 5:-1:1
        for j = i-1:-1:1
            disp(pVals(pr, j, i));
        end
    end
    disp("~~~")
end

%% 
% disp('v1')
% [P,ANOVATAB] = kruskalwallis(squeeze(allBootDists(1, :, :))')
% disp('ppa')
% [P,ANOVATAB] = kruskalwallis(squeeze(allBootDists(2, :, :))')

%%