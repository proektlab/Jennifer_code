%% function: better phase probablility plot generator
%gives MI, a statistic that shows how entrained a neuron is to the deltawave,
%phase with max amplitude, can also plot phase vs probability
%input: delta data - lfp filtered for delta, times x channel x trial
%       trialSpiketrain - cell matrix that is 1 x number of trials
%       cointaining timepoints of spikes in each trial. spikes should have
%       time relative to stimulus.
%       figureON - true for figure, false for no figure
%       limits (optional) - timepoints (relative to trial start) to look at
%       example: [100 900] for 100 to 900ms. This function assumes
%       stimulus is at 1000ms
%       nbin (optional)- number of bins, default 20
%       deltaStartTime (optional) time the delta trial starts, default
%       1000ms
%output:  trialMeanSpikes, mean number of spikes in each phase bin
%         maxAmp - bin with most spikes
%         MI - modulation index, measures how uniform it is
%jl 6/17/21 copied, 7/27/21 edited
function [trialMeanSpikes, maxAmp, MI] = betterPhaseProb(deltaData, trialSpiketrain, figureOn, limits, nbin, deltaStartTime)

if nargin < 6
    deltaStartTime = 1000; %milliseconds before stimulus
end
if nargin < 5
    nbin = 20; % number of phase bins
end
if nargin < 4
    limits = [1000, 1900]; %milliseconds relative to trial start
end
numTrials = numel(trialSpiketrain);

ch = 1;
delSig = permute(deltaData, [2, 3, 1]);
%find phase of the delta filtered data
for tr = 1:numTrials
    phaseDelta(ch,tr,:) =  angle(hilbert(squeeze(delSig(ch,tr,:))));
end

%make bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi+(j-1)*winsize;
end

trialSpikes = zeros(numTrials, nbin);
%for each trial
for trial = 1:numTrials
    %for each bin
    for bin = 1:nbin
        %find number of spikes in  trialSpiketrain that are in this bin
        binLocs = find(phaseDelta(1, trial, :) <  position(bin)+winsize &...
            phaseDelta(1, trial, :) >=  position(bin));
        binLocs = binLocs(binLocs > (limits(1)) & binLocs < (limits(2))) - deltaStartTime; %adjustment with trialspiketime
        trialSpikes(trial, bin) = numel(intersect(binLocs, floor(trialSpiketrain{trial}))); %number of spikes that are also in these phases
    end
end

%average over trials
trialMeanSpikes = mean(trialSpikes, 1);

%generate bin with largest amplitude
[~, ind] = max(trialMeanSpikes);
maxAmp = position(ind);

%calculate MI
MI=(log(nbin)-(-nansum((trialMeanSpikes/nansum(trialMeanSpikes)).*log((trialMeanSpikes/nansum(trialMeanSpikes))))))/log(nbin);

if(figureOn)
    bar(position, trialMeanSpikes/sum(trialMeanSpikes));
    %labels
    xlabel("Phase")
    ylabel("Actual Spike Probability")
end
end