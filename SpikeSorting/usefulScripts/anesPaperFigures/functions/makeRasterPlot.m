%function makeRasterPlot(trialSpiketrain, deltaPlot, gammaPlot, histogramOn, color)
%this function plots a raster plot, with gamma and delta overlays optional
%inputs:    trialSpiketrain - 1xnumTrials cell vector containing spiketimes of
%           every trial
%           deltaPlot - values to plot for a delta overlay. empty if no
%           plot.
%           gamma plot - values for a gamma overlay, empty if no plot
%           histogramOn - whether or not there'll be a psth
%           color - color of plots
%           stimTime - time of stimulus relative to start of trial
%jl 8/5/21
function makeRasterPlot(trialSpiketrain, deltaPlot, gammaPlot, histogramOn, color, stimTime)
if nargin < 6
    stimTime = 1000;
end
if nargin < 5
    color = 'b';
end

if nargin < 4
    histogramOn = 0;
end
if nargin < 3
    gammaPlot = [];
end
if nargin < 2
    deltaPlot = [];
end

x = [];
y = [];
%temp
trialsToUse = 1:numel(trialSpiketrain);
sizeDelta = size(deltaPlot, 1);
plotEndpoint = sizeDelta - 1 - stimTime;

%get x and y values
for trial = trialsToUse
    x = [x, trialSpiketrain{trial}];
    y = [y, trial * ones(size(trialSpiketrain{trial}))];
end

barPlotNums = zeros(1, 100);
%histogram bin numbers
for bin = 1:100
    currBinLims = [-1000 + 30*(bin - 1), -970 + 30*(bin - 1)];
    barPlotNums(bin) = nnz(x < currBinLims(2) & x > currBinLims(1));
end
%delta overlay
if (~isempty(deltaPlot))
    if(size(deltaPlot, 2) > 1)
        numDeltas = size(deltaPlot, 2);
        for iPlot = 1:floor(numDeltas/10):numDeltas
            plot(-stimTime:plotEndpoint, deltaPlot(:, iPlot)./max(deltaPlot(:, iPlot))*5 + 20 + iPlot, 'color', [0.8500 0.3250 0.0980]);
            hold on
        end
    else
        plot(-stimTime:plotEndpoint, deltaPlot./max(deltaPlot)*20 + 25, 'color', [0.8500 0.3250 0.0980]);
        hold on
    end
end
%gamma overlay
if(~isempty(gammaPlot))
    plot(-stimTime:plotEndpoint, averageGamma/5 - 25, 'color', [0.4660 0.6740 0.1880]);
    plot(-stimTime:plotEndpoint, gammaAmp/5 - 25, 'color', [0.4660 0.6740 0.1880], "LineWidth", 1);
end

scatter(x, y + 20, 15, color , 'filled', 's', 'MarkerFaceAlpha',.75);
ylim([-10,size(trialSpiketrain,2)+20]);
yticks(20:20:120);
yticklabels({0:20:100});
xlim([-500 1000]);
hold on;

xline(0);

%psth
if(histogramOn)
    bar(-970:30:2000, barPlotNums/max(barPlotNums) * 20, 'FaceAlpha', .5, 'FaceColor', color)
end
hold off

xlabel("Time (ms)")

end