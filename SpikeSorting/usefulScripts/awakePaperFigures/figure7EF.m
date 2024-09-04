%generates phase/probabilty images for all entrained awake cells, sorted by
%layer

%% load things
load('poissonSurrogateDataAwakeLED.mat')

dataPath = 'Z:\adeeti\JenniferHelen\SpikeSortingResults\LEDAwakeOnly\';
experimentFolders = dir([dataPath, 'Output*']);

layerCats = {1, [2, 3], 4, 5, 6, [9, 10, 11, 12], [13, 14]};
%% get some info
%start a matrix for each layer category
allPhaseProbs = cell(1,7);
allMaxAmps = cell(1,7);
someChosenSpikes = [];
%for each experiment
for i = 1:numel(experimentFolders)
    %load post kilo sort file
    cd([dataPath, experimentFolders(i).name, '\Data']); %go to that data folder
    load(['DataFromPhy-', experimentFolders(i).name(8:end), '.mat']) %load  datafiles
    load('usefulVars')
    %possible locations of delta
    depthName = ['Z:\adeeti\ecog\iso_awake_VEPs\depthMice\', experimentFolders(i).name(8:11)];
    goodName = ['Z:\adeeti\ecog\iso_awake_VEPs\goodMice\', experimentFolders(i).name(8:11)];
    %go to correct layer data location
    if exist(depthName, 'file')
        load([depthName '\CSD_FiltData\', experimentFolders(i).name(8:end), '.mat']...
            , 'csdPr1_filtSigDelta', 'csdPr2_filtSigDelta') %load probe delta
    else
        if exist(goodName, 'file')
            load([goodName '\CSD_FiltData\', experimentFolders(i).name(8:end), '.mat']...
                , 'csdPr1_filtSigDelta', 'csdPr2_filtSigDelta') %load probe delta
        end
    end
    csdDelta = horzcat(csdPr1_filtSigDelta, csdPr2_filtSigDelta);
    numUnits = numel(channel);
    %for each unit
    for unit = 1:numUnits
        %if it's significant
        usableUnit = ~isnan(significance{i}(unit))&&significance{i}(unit);
        if(usableUnit)
            %calculate phase prob
            currChan =  trueLoc(unit)-64;
            while(nnz(isnan(csdDelta(:, currChan, :))))
                currChan = currChan - 1;
            end
            currLayerCat = find(cell2mat(cellfun(@(x) any(x == unitLayers(unit)), layerCats, 'UniformOutput', false)));
            [currTrialMeanAmp, maxAmp] = betterPhaseProb(csdDelta(:, currChan, :), trialSpiketrain(unit,2:end), 0);
            if(~isempty(currLayerCat))
                allPhaseProbs{currLayerCat}(end + 1,:) = currTrialMeanAmp/sum(currTrialMeanAmp); %and add it to its relevant layer
                allMaxAmps{currLayerCat}(end +1) = maxAmp;
            else
                disp('nyeh? (csd delta has nans?)')
            end
            
            %determine which ones I plotted
            if currLayerCat == 5
                someChosenSpikes(end+1) = 0;
                if unit == 25 && i == 3
                    someChosenSpikes(end) = 1;
                end
                if unit == 27 && i == 3
                    someChosenSpikes(end) = 2;
                end
            end
        end
    end
end

%% units to use
%CB34 unit 17 phase 1.26
%CB34 unit 12 phase -2.51
%both are V1
% figure
% counter = 1;
% for i = 1:numUnits
%     if unitLayerCat(i) > 5
%         continue
%     end
%     currChan =  trueLoc(i)-64;
%     while(nnz(isnan(csdDelta(:, currChan, :))))
%         currChan = currChan - 1;
%     end
%     [currTrialMeanAmp, maxAmp] = betterPhaseProb(csdDelta(:, currChan, :), trialSpiketrain(i,2:end), 0);
%     subplot(5, 5, counter)
%     counter = counter + 1;
%     makeRasterPlot(trialSpiketrain(i, 2:end));
%     title([num2str(maxAmp),' ', num2str(unitLayers(i)),' ', num2str(i)]);
% end
% sgtitle("CB26")

%% sort some things
for layerCat = 2:7
    [~, sortInd] = sort(allMaxAmps{layerCat});
    sortedPhaseProbs{layerCat} = allPhaseProbs{layerCat}(sortInd, :);
end

%for finding the ones to pick out
[~, sortInd] = sort(allMaxAmps{5});
newLocs = someChosenSpikes(sortInd);
%% make some figures
%get some colormaps that are pretty
originalColors = [0 0 1];
v1Colormap = zeros(256, 3);
for col = 1:3
    v1Colormap(:, col) = linspace(1, originalColors(col), 256);
end

originalColors = [1 0 0];
ppaColormap = zeros(256, 3);
for col = 1:3
    ppaColormap(:, col) = linspace(1, originalColors(col), 256);
end

layerNames = ["", "V1 Layer 2/3", "V1 Layer 4", "V1 Layer 5", "V1 Layer 6", "Upper PPA", "Lower PPA"];
figure
%for each layer category
for layerCat = 2:5
    %subplot
    subplot(1, 4, layerCat - 1);
    %image the phase probs
    imagesc(sortedPhaseProbs{layerCat});
    colormap(v1Colormap)
    title(layerNames(layerCat))
    xlabel("phase")
    ylabel("unit")
    xticks(0:5:20)
    xticklabels({'-pi', '-pi/2', '0', 'pi/2', 'pi'});
    colorbar
end

sgtitle("Phase Probability of All Significant V1 Units")

figure
%for each layer category
for layerCat = 6:7
    %subplot
    subplot(1, 2, layerCat - 5);
    %image the phase probs
    imagesc(sortedPhaseProbs{layerCat});
    colormap(v1Colormap)
    title(layerNames(layerCat))
    xlabel("phase")
    ylabel("unit")
    xticks(0:5:20)
    xticklabels({'-pi', '-pi/2', '0', 'pi/2', 'pi'});
end

sgtitle("Phase Probability of All Significant PPA Units")

%% stack some things
allV1s = vertcat(sortedPhaseProbs{2}, sortedPhaseProbs{3}, sortedPhaseProbs{4}, sortedPhaseProbs{5});
v1Borders = cumsum(layerSig);
v1Borders = v1Borders(2:4);

allPPAs = vertcat(sortedPhaseProbs{6}, sortedPhaseProbs{7});
ppaBorder = layerSig(6);

ff = figure('Color', 'w', 'Renderer', 'Painters');

%subplot(1, 2, 1)
ax1 = axes('Position',[0.1 0.1 0.35 0.8]);
imagesc(ax1, allV1s);
hold on
for i = 1:3
    yline(v1Borders(i) + .5, 'LineWidth', 1);
end
%colorbar
colormap(ax1, ppaColormap);
c = colorbar;
c.Label.String = 'Fraction of Spikes';
caxis([0.03, .15]); %colorbar limits
%labels
title('V1')
xlabel("Phase of Delta")
ylabel("Unit ID")
xticks(0:5:20)
xticklabels({'-pi', '-pi/2', '0', 'pi/2', 'pi'});

%subplot(1, 2, 2)
ax2 = axes('Position',[0.6 0.1 0.35 0.8]);
imagesc(ax2, allPPAs);
hold on
for i = 1
    yline(ppaBorder(i) + .5, 'LineWidth', 1);
end
colormap(ax2, v1Colormap);
c = colorbar;
c.Label.String = 'Fraction of Spikes';
caxis([0.03, .12]); %colorbar limits
%labels
title('PPA')
xlabel("Phase of Delta")
ylabel("Unit ID")
xticks(0:5:20)
xticklabels({'-pi', '-pi/2', '0', 'pi/2', 'pi'});

%% 
figure('color', 'white')
rectLocs = cellfun(@numel, sortedPhaseProbs)/20;
rectLocs = sum(rectLocs(1:4));
rectLocs = rectLocs + find(newLocs);

imagesc(allV1s);
hold on
for i = 1:3
    yline(v1Borders(i) + .5, 'LineWidth', 1);
end
rectangle('Position', [0, (rectLocs(1)- .5), 21, 1], 'EdgeColor', 'b', 'LineWidth', 1) 

rectangle('Position', [0, (rectLocs(2)- .5), 21, 1], 'EdgeColor', 'b', 'LineWidth', 1) 

%colorbar
colormap(ppaColormap);
c = colorbar;
c.Label.String = 'Fraction of Spikes';
caxis([0.03, .15]); %colorbar limits
%labels
title('V1')
xlabel("Phase of Delta")
ylabel("Unit ID")
xticks(0:5:20)
xticklabels({'-pi', '-pi/2', '0', 'pi/2', 'pi'});
%% to excel
figureE = allV1s;
figureE = array2table(figureE);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7EPhaseProbV1.xlsx';
writetable(figureE,filename,'Sheet',1,'Range','A1')

figureF = allPPAs;
figureF = array2table(figureF);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7FPhaseProbPPA.xlsx';
writetable(figureF,filename,'Sheet',1,'Range','A1')