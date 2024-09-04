% function plotWaveformsFilt(mouseName, chosenSpike, color)
%finds the waveforms of a given unit by filtering its wideband data and
%plots them
%inputs:    mouseName - name of trial ('CB34_Awake_LED100")
%           expType - type of experiment ('Awake', 'LowIso', or 'Ket')
%           chosenSpike - index of spike
%           color - color of plot
%           createPlot - whether the plot should be created
% jl 10/12/23 

function actualWaveforms = plotWaveformsFilt(mouseName, expType, chosenSpike, color, createPlot)
if nargin < 4
    color = 'blue';
end

if nargin < 5
    createPlot = true;
end

shortName = mouseName(1:4);

%load appropriate file
[~, trueLoc, ~, spiketrain] = getUsefulVars_allData(mouseName, expType);

% filter
myExperiment = ['Z:\code\Jennifer_code\SpikeSorting\usefulScripts\anesPaperFigures\allData\Representative_Neurons_pl2\',...
    mouseName, '.pl2'];

[MUAData, ~, ~, ~, ~, ~] = ExtractPL2_MUA_from_WB_Plexon(myExperiment,...
    [], 1, 325, 5000, 40000, trueLoc(chosenSpike));

% get the *actual* waveforms
chosenSpikeChan = trueLoc(chosenSpike); %%find spike location
numSpikes = numel(spiketrain{chosenSpike});
actualWaveforms = zeros(161, numSpikes);
for spike = 1:numSpikes
    spikeInd = [ceil(spiketrain{chosenSpike}(spike)/1000*40000 - 80), ceil(spiketrain{chosenSpike}(spike)/1000*40000 + 80)];
    actualWaveforms(:, spike) = MUAData(:, spikeInd(1):spikeInd(2));
end

if createPlot
    % plot
    for i = 1:10:numSpikes
        x = -2:1/40:2;
        y = detrend(actualWaveforms(:, i)');
        plot(x, y, 'Color', color)
        hold on
    end
    plot(x, mean(actualWaveforms, 2), 'k');
    xlabel('Time (ms)')
    ylabel('Voltage');
end
end