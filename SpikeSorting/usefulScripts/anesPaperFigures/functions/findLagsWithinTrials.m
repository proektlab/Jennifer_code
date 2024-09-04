%function allLags = findLagsWithinTrials(trialSpiketrain1, trialSpiketrain2, timeLimit) 
%finds all the lags between spiketrain1 and spiketrain2 within a trial
%output: all time lags between trialspiketrain1 and trialspiketrain2
%input
%   trialSpiketrain1 - spiketrain, cell vector in 1 x number of trials,
%   with spike times relative to stimulus in every trial. 
%   trialSpiketrain2 - format as above, will have its spike times
%   subtracted from the spiketimes of 1
%   timeLimit - 2 numbers containing lower limit and upper limit of spike
%   times, such as [-1000 0]
function allLags = findLagsWithinTrials(trialSpiketrain1, trialSpiketrain2, timeLimit) 
    if nargin == 2
        timeLimit = [-1000, 2000];
    end
    allLags = [];
    numTrials = numel(trialSpiketrain1);
    %for each trial
    for trial = 1:numTrials
        currST1 = trialSpiketrain1{trial}; 
        currST2 = trialSpiketrain2{trial}; 
        limitedST1 = currST1(currST1<timeLimit(2)&currST1>timeLimit(1)); 
        limitedST2 = currST2(currST2<timeLimit(2)&currST2>timeLimit(1));
        for spikeTime = 1:numel(limitedST2)
            %find the difference between it and every spike in chosen2
            %and add it to prestimlags
            allLags = [allLags, limitedST1 - limitedST2(spikeTime)];
        end
    end
end