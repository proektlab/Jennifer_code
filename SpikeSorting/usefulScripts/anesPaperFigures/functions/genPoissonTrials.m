%geneates some poisson trials from a spike probability vector
%inputs:    probs -  vector of spike probabilities for each timepoint
%           numTrials - number of trials to generate
%           surrogate trials - cell vector of spiketrains for each poisson
%           trial
%           timepoints - timepoints relative to stimulus over which to
%           generate trials
%outputs: surrogateTrials - 1xnumTrials cell containing spike times of
%         surrogate neuron for each trial
function surrogateTrials = genPoissonTrials(probs, numTrials, timepoints)
surrogateTrials = cell(1,numTrials);
%for each trial
for n = 1:numTrials
    %for each timestamp
    for m = 1:numel(probs)
        %generate a random number
        %if it's above corresponding probability
        if(rand < probs(m))
            %add a spike there
            surrogateTrials{n} = [surrogateTrials{n}, timepoints(m)];
        end
    end
end
end