
%This takes a representative iso, ket and awake neuron, and time warps each
%trial such that the phases of the delta wave in the trial are aligned.
%Then, it plots the time-warped raster plot
%% make a group of figures
ff = figure('Color', 'w', 'Renderer', 'Painters');
chosenMice = [27,  37,  40];
chosenUnit = [25,  7,  13];
chosenType = {'Awake', 'LowIso', 'Ket100'};
anesColors = {[0.6350 0.0780 0.1840], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
smooshMean = zeros(2, 3); %for recording the smoosh

for i = 1:3 %for each condition
    %load relevant files
    %get the csd 
     currMouseName = ['CB', num2str(chosenMice(i)), '_', chosenType{i}, '_LED100'];
    if strcmp(chosenType{i}, 'Ket')
        currMouseName = ['CB', num2str(chosenMice(i)), '_', chosenType{i}, '100_LED100'];
    end
    [trialSpiketrain, trueLoc, unitLayercat, ~] = getUsefulVars_allData(currMouseName, chosenType{i});
    csdDelta = getCSDDelta_allData(currMouseName, chosenType{i});

    %subplot
    subplot(3, 3, i)
    %normal raster
    makeRasterPlot(trialSpiketrain(chosenUnit(i), :), squeeze(csdDelta(:, trueLoc(chosenUnit(i))-64, :)), [], true, anesColors{i})
    title([chosenType(i), 'Raster']);
    
    %subplot
    subplot(3, 3, 3 + i)
    %pre stim squished raster
    title('Prestim Time Warped Raster')
    preSmoosh = makeSquishedPhaseRaster(trialSpiketrain(chosenUnit(i), :), squeeze(csdDelta(:, trueLoc(chosenUnit(i))-64, :)), [], true, anesColors{i}, 100:900);
    smooshMean(1, i) = mean(preSmoosh);
    
    %subplot
    subplot(3, 3,6 + i)
    %poststim squished raster
    title('Poststim Time Warped Raster')
    postSmoosh = makeSquishedPhaseRaster(trialSpiketrain(chosenUnit(i), :), squeeze(csdDelta(:, trueLoc(chosenUnit(i))-64, :)), [], true, anesColors{i}, 1100:2000);
    smooshMean(2, i) = mean(postSmoosh);
end

%add some labels 
subplot(3, 3, 1);
xlabel("Time (ms)")
ylabel("Trial");

subplot(3, 3, 2);
xlabel("Time (ms)")

subplot(3, 3, 3);
xlabel("Time (ms)")

subplot(3, 3, 4);
ylabel("Trial");

subplot(3, 3, 7);
ylabel("Trial");
xlabel("Time (Adjusted)");

subplot(3, 3, 8);
xlabel("Time (Adjusted)")

subplot(3, 3,9);
xlabel("Time (Adjusted)")