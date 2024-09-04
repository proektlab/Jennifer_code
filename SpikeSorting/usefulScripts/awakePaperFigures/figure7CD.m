%%yet another thing that makes rasters
% but this one uses trial spiketrain, now that it's been made
%(generates rasters for awake paper)
%(can use makeRasterPlot.m function instead)
%6/30/21 jl

%% chosen dataset: CB34, unit 19 for v1
mouseName = 'CB34';
trialName = 'CB34_Awake_LED100';


%load appropriate file
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\DataFromPhy-', trialName])
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\usefulVars'])
%load filtered stuff
load(['Z:\adeeti\ecog\iso_awake_VEPs\depthMice\', mouseName, '\CSD_FiltData\', trialName, '.mat'], 'csdPr1_filtSigGamma', 'csdPr1_filtSigDelta')

%%
figure
chosenUnit = 19;
x = [];
y = [];
for trial = 1:100
    x = [x, trialSpiketrain{chosenUnit, trial}];
    y = [y, trial * ones(size(trialSpiketrain{chosenUnit, trial}))];
end

barPlotNums = zeros(1, 100);
for bin = 1:100
    currBinLims = [-1000 + 30*(bin - 1), -970 + 30*(bin - 1)];
    barPlotNums(bin) = nnz(x < currBinLims(2) & x > currBinLims(1));
end

averageDelta = squeeze(mean(csdPr1_filtSigDelta(:, trueLoc(chosenUnit) - 64, :), 3));
averageGamma = squeeze(mean(csdPr1_filtSigGamma(:, trueLoc(chosenUnit) - 64, :), 3));
gammaAmp = abs(hilbert(averageGamma));

scatter(x, y + 20, 25,[0 0.4470 0.7410] , 'filled', 's', 'MarkerFaceAlpha',.25);
hold on;
plot(-1000:2000, averageDelta/5 + 25, 'color', [0.8500 0.3250 0.0980]);
plot(-1000:2000, averageGamma/5 - 25, 'color', [0.4660 0.6740 0.1880]);
plot(-1000:2000, gammaAmp/5 - 25, 'color', [0.4660 0.6740 0.1880], "LineWidth", 1);
bar(-970:30:2000, barPlotNums/max(barPlotNums) * 20, 'FaceAlpha', .5, 'FaceColor', [0 0.4470 0.7410] )

xlabel("Time (ms)")
ylim([-60, 120])
title("V1 Unit Raster Plot")

%% to excel
%bar
figureC = [-970:30:2000; barPlotNums];
figureC = array2table(figureC);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7CBar.xlsx';
writetable(figureC,filename,'Sheet',1,'Range','A1')

%% chosen unit: ppa

mouseName = 'CB37';
trialName = 'CB37_Awake_LED100';


%load appropriate file
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\DataFromPhy-', trialName])
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\usefulVars'])
%load filtered stuff
%load(['Z:\adeeti\ecog\iso_awake_VEPs\depthMice\', mouseName, '\prLFP_FiltData\', trialName, '.mat'], 'pr1_filtSigDelta', 'pr2_filtSigDelta')
load(['Z:\adeeti\ecog\iso_awake_VEPs\depthMice\', mouseName, '\CSD_FiltData\', trialName, '.mat'], 'csdPr2_filtSigGamma', 'csdPr2_filtSigDelta')
%%
%subplot(1, 2, 2)
ff = figure('Color', 'w', 'Renderer', 'Painters');
chosenUnit = 3;
x = [];
y = [];
for trial = 1:100
    x = [x, trialSpiketrain{chosenUnit, trial}];
    y = [y, trial * ones(size(trialSpiketrain{chosenUnit, trial}))];
end

barPlotNums = zeros(1, 100);
for bin = 1:100
    currBinLims = [-1000 + 30*(bin - 1), -970 + 30*(bin - 1)];
    barPlotNums(bin) = nnz(x < currBinLims(2) & x > currBinLims(1));
end

averageDelta = squeeze(mean(csdPr2_filtSigDelta(:, trueLoc(chosenUnit) - 96, :), 3));
averageGamma = squeeze(mean(csdPr2_filtSigGamma(:, trueLoc(chosenUnit) - 96, :), 3));
gammaAmp = abs(hilbert(averageGamma));

scatter(x, y + 20, 25,'b' , 'filled', 's', 'MarkerFaceAlpha',.25);
hold on;
plot(-1000:2000, averageDelta/5 + 25, 'color', [0.8500 0.3250 0.0980]);
plot(-1000:2000, averageGamma/5 - 25, 'color', [0.4660 0.6740 0.1880]);
plot(-1000:2000, gammaAmp/5 - 25, 'color', [0.4660 0.6740 0.1880], "LineWidth", 1);
bar(-970:30:2000, barPlotNums/max(barPlotNums) * 20, 'FaceAlpha', .5, 'FaceColor', 'b')
xline(0, "--g")

xlabel("Time (ms)")
ylim([-60, 120])
title("PPA Unit Raster Plot")

%% to excel
%bar
figureD = [-970:30:2000; barPlotNums];
figureD = array2table(figureD);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7DBar.xlsx';
writetable(figureD,filename,'Sheet',1,'Range','A1')

%% show the spike waveforms (broken :( )
% chosenUnit = 21;
% figure()
% for i = 1:10%:size(waveforms{chosenUnit}, 2)
%     x = 1:103;
%     y = detrend(waveforms{chosenUnit}(:, i)');
%     plot(x, y, 'b')
%     hold on
%     
% end
% plot(averageWF{chosenUnit}, 'r');
% ylim([-5, 5])
% xlabel('time')
% ylabel('voltage')
%% testing some things
% 
% figure()
% [sorted, index] = sort(trueLoc);
% sorted = sorted(end-6:end);
% index = index(end-6:end);
% for i = 1:7
%     h(i)=subplot(8, 1, i)
%     for j = 1:100
%         plot(-1000:2000, csdPr2_filtSigDelta(:, trueLoc(index(i))-96, j), 'color', [0 0.4470 0.7410]);
%         hold on
%     end
% end
% linkaxes(h, 'xy');