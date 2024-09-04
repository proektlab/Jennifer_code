%This displays some waveforms and rasters of a representative neuron from
%V1 and PPA, and under each condition


%chosen mice
anesColors = {[0.6350 0.0780 0.1840], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
chosenUnit = [20, 19, 23];
chosenMice = [40, 34, 37];
chosenType = {'Awake', 'LowIso', 'Ket100'};
chosenLoc = ["PPA", "PPA", "PPA"];

%figure
waveformFig = figure('Color', 'w', 'Renderer', 'Painters');

%figure
rasterFig = figure('Color', 'w', 'Renderer', 'Painters');

%loop
%for each representative neuron
for ind = 1:numel(chosenUnit)
    currColor = anesColors{ind}; %each condition is color coded
    currMouseName = ['CB', num2str(chosenMice(ind)), '_', chosenType{ind}, '_LED100'];
    if strcmp(chosenType{ind}, 'Ket')
        currMouseName = ['CB', num2str(chosenMice(ind)), '_', chosenType{ind}, '100_LED100'];
    end
    
    figure(waveformFig)
    subplot(2, 3, ind)
    %plot the waveform
    plotWaveformsFilt(currMouseName, chosenType{ind}, chosenUnit(ind), currColor)
    ax = gca;
    %title
    title([currMouseName, ' Unit ' num2str(chosenUnit(ind)),...
        ', ', convertStringsToChars(chosenLoc(ind)), ' Neuron'],'Interpreter','none');
    
    %figure(rasterFig)
    subplot(2, 3,ind+3)
    %load spiketrain data
    [trialSpiketrain, trueLoc, ~, ~] = getUsefulVars_allData(currMouseName, chosenType{ind});
    %load csd data
    csdDelta = getCSDDelta_allData(currMouseName, chosenType{ind});
    %plot figure
    makeRasterPlot(trialSpiketrain(chosenUnit(ind), :), squeeze(mean(csdDelta(:, trueLoc(chosenUnit(ind)) - 64,:), 3)), [], true, currColor);
    ax = gca;
    %title
    title([currMouseName, ' Unit ' num2str(chosenUnit(ind)),...
        ', ', convertStringsToChars(chosenLoc(ind)), ' Raster Plot'],'Interpreter','none');
end
