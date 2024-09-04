%just trying to get some spike waveforms ;_;
%(generate spike waveforms for awake paper)
%(can use plotWaveformsFilt.m function instead)

%% filter

myExperiment = 'D:\adeeti\Data4KiloSort\LEDFlashes\CB34\CB34_Awake_LED100.pl2';

[MUAData, ogSampRate, eventTimes, fullTraceTime, plexInfoStuffs, fragInfo] =...
    ExtractPL2_MUA_from_WB_Plexon(myExperiment);

%% grab spike times

mouseName = 'CB34';
trialName = 'CB34_Awake_LED100';


%load appropriate file
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\DataFromPhy-', trialName])
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\usefulVars'])

%% get the *actual* waveforms
chosenSpike = 19;
chosenSpikeChan = trueLoc(chosenSpike); %hopefully this is the right one?
numSpikes = numel(spiketrain{chosenSpike});
actualWaveforms = zeros(161, numSpikes);
for spike = 1:numSpikes
    spikeInd = [ceil(spiketrain{chosenSpike}(spike)/1000*40000 - 80), ceil(spiketrain{chosenSpike}(spike)/1000*40000 + 80)];
    actualWaveforms(:, spike) = MUAData(chosenSpikeChan - 64, spikeInd(1):spikeInd(2));
end
% 
%% plot
figure('Color', 'w', 'Renderer', 'Painters')
%subplot(1, 2, 1);
for i = 1:10:numSpikes
    x = -2:1/40:2;
    y = detrend(actualWaveforms(:, i)');
    plot(x, y, 'Color', [0 0.4470 0.7410])
    hold on
end
plot(x, mean(actualWaveforms, 2), 'k');
xlabel('Time (ms)')
ylabel('Voltage');
title('Example of V1 Neuron')

%% save into excel file
figureA = [x; mean(actualWaveforms, 2)'];
figureA = array2table(figureA);

 filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7AWaveform.xlsx';
 writetable(figureA,filename,'Sheet',1,'Range','A1')
% %% do ksdensity for v1
% figure(2)
% tiledlayout(1, 2);
% ax1 = nexttile;
% clear probDens
%  for ch = 1%:numSpikes
%     useData = actualWaveforms;
%     maxLim = .25%max(useData(:));
%     if isnan(maxLim)
%         maxLim =0;
%     elseif maxLim>750
%         maxLim = 750;
%     end
%     minLim = -.30; %min(useData(:));
%     if isnan(minLim)
%         minLim =0;
%     elseif minLim<-750
%         minLim = -750;
%     end
%     pts = linspace(minLim, maxLim, 500);
%         if sum(isnan(useData(:)))==numel(useData)
%             probDens(ch,:,:) = useData(1:50,:);
%         else
%             for t = 1:161
%                 [f,xi] =  ksdensity(useData(t, :),pts);
%                 probDens(:, t) = f;
%             end
%         end
%  end
%  imagesc(flip(probDens));
%  
%  %make a neat looking color map
%  originalColors = [0 0.4470 0.7410];
%  v1Colormap = zeros(256, 3);
%  for col = 1:3
%      v1Colormap(:, col) = linspace(1, originalColors(col), 256);
%  end
%  
%  colormap(ax1, v1Colormap);
%  
%  %labels
%  yticklabels(flip(round(xi(1:50:500), 2)))
%  ylabel("Voltage")
%  xlabel("Time (ms)")
%% filter

myExperiment = 'D:\adeeti\Data4KiloSort\LEDFlashes\CB37\CB37_Awake_LED100.pl2';

[MUAData2, ogSampRate2, eventTimes2, fullTraceTime2, plexInfoStuffs2, fragInfo2] =...
    ExtractPL2_MUA_from_WB_Plexon(myExperiment, [], 32, 325, 5000, 40000, 96);

%% grab spike times

mouseName = 'CB37';
trialName = 'CB37_Awake_LED100';


%load appropriate file
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\DataFromPhy-', trialName])
load(['D:\adeeti\KiloSort_results\LEDAwakeOnly\Output_', trialName, '\Data\usefulVars'])

%% get the *actual* waveforms
chosenSpike = 3;
chosenSpikeChan = trueLoc(chosenSpike); %hopefully this is the right one?
numSpikes = numel(spiketrain{chosenSpike});
actualWaveforms = zeros(161, numSpikes);
for spike = 1:numSpikes
    spikeInd = [ceil(spiketrain{chosenSpike}(spike)/1000*40000 - 80), ceil(spiketrain{chosenSpike}(spike)/1000*40000 + 80)];
    actualWaveforms(:, spike) = MUAData2(chosenSpikeChan - 96, spikeInd(1):spikeInd(2));
end

%actualWaveforms = waveforms{chosenSpike};

%% plot
ff = figure('Color', 'w', 'Renderer', 'Painters');
%subplot(1, 2, 2);
for i = 1:10:numSpikes
    x = -2:1/40:2;
    y = detrend(actualWaveforms(:, i)');
    plot(x, y, 'Color', 'b')
    hold on
end
plot(x, mean(actualWaveforms, 2), 'k');
xlabel('Time (ms)')
ylabel('Voltage');
title('Example of PPA Neuron')

%% save into excel file
figureB = [x; mean(actualWaveforms, 2)'];
figureB = array2table(figureB);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figure7BWaveform.xlsx';
writetable(figureB,filename,'Sheet',1,'Range','A1')

%% do ksdensity for ppa
% figure(2)
% ax2 = nexttile;
% clear probDens
%  for ch = 1%:numSpikes
%     useData = actualWaveforms;
%     maxLim = .25;%max(useData(:));
%     if isnan(maxLim)
%         maxLim =0;
%     elseif maxLim>750colo
%         maxLim = 750;
%     end
%     minLim = -.25; %min(useData(:));
%     if isnan(minLim)
%         minLim =0;
%     elseif minLim<-750
%         minLim = -750;
%     end
%     pts = linspace(minLim, maxLim, 500);
%         if sum(isnan(useData(:)))==numel(useData)
%             probDens(ch,:,:) = useData(1:50,:);
%         else
%             for t = 1:161
%                 [f,xi] =  ksdensity(useData(t, :),pts);
%                 probDens(:, t) = f;
%             end
%         end
%  end
%  imagesc(flip(probDens));
%  
%  %make a neat looking color map
%  originalColors = [0.6350 0.0780 0.1840];
%  ppaColormap = zeros(256, 3);
%  for col = 1:3
%      ppaColormap(:, col) = linspace(1, originalColors(col), 256);
%  end
%  
%  colormap(ax2, ppaColormap);
%  
%  %labels
% yticklabels(flip(round(xi(1:50:500), 2)))
%  ylabel("Voltage")
%  xlabel("Time (ms)")