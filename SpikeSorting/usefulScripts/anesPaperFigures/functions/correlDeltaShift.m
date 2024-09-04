%this function finds the correlation of two trials of delta, one which has
%been shifted and compressed
%inputs
%       push: adds some padding to the beginning of delta1
%       smoosh: factor by which the data is timeshifted
%       delta1: first trace of delta
%       delta2: second trace of delta, which will be transformed
%       timePoints_a: timepoints of delta1 to use for analysis
%       timepoints_b: timepoitns of delta2 to use for analysis
%       trialSpiketrain: (optional) spiketrain to transform using the push
%       and smoosh values
%outputs
%       correls: correlation between delta1 and the transformed delta2
%       newSpiketimes: spiketrain transformed using the push and smoosh
%       values
%       corr1: delta1 at the selected timepoints
%       corr2: delta2 at the selected timepoints transformed with the push
%       and smoosh values
%       timepoints2: timepoints for corr2

function [correls, newSpiketimes, corr1, corr2, timepoints2] = correlDeltaShift(push, smoosh, delta1, delta2, timePoints_a, timePoints_b, trialSpiketrain1)
    if nargin < 7
       trialSpiketrain1 = [];
    end
   
    
    timepoints1 = timePoints_a;
    timepoints2 = linspace(timePoints_b(1), timePoints_b(end), numel(timePoints_a));
    timepoints2 = (timepoints2 - timePoints_b(1))*smoosh +timePoints_b(1); %scale by smoosh value
    timepoints2 = timepoints2+ push; %add push value to get new timepoints
    
    %get delta at those timepoints, in the case of corr2 perform
    %interpolation
    corr1 = delta1(timepoints1);
    corr2 = interp1(1:numel(delta2), delta2, timepoints2); 
    
    %get correlation value
    if numel(corr1) ~= numel(corr2)
        disp('???')
    end
    correls = corr(corr1, corr2');  
    
    if isempty(trialSpiketrain1) 
        newSpiketimes = [];
        return
    end
    
    
    %shift spiketimes so that they can be plotted with the delta
    newSpiketimes = trialSpiketrain1(trialSpiketrain1 < timepoints2(end)  &  trialSpiketrain1 > timepoints2(1));
    newSpiketimes = (newSpiketimes - push - timePoints_b(1) )/smoosh + timePoints_b(1);
    
    %scale spiketimes properly
    newSpiketimes = (newSpiketimes - timePoints_b(1))/numel(timePoints_b)*numel(timePoints_a);
    

    
end

