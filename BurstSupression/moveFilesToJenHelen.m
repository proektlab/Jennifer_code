dirIn = 'Z:\adeeti\ecog\matIsoPropMultiStim\';
dirOut = 'Z:\adeeti\JenniferHelen\';

cd(dirIn)

flashIndex = [0 inf];
whiskIndex = [inf 0];

%[myFavoriteExp] = findMyExpMulti(dataMatrixFlashes, exp, drugType, conc, stimIndex,  numStim, typeTrial, forkPos)

%% high iso flashes
[myFavoriteExp] = findMyExpMulti(dataMatrixFlashes, [], 'iso', 1.2, flashIndex);
dirOutIsoFlashes = [dirOut, 'Iso_flashes\'];

mkdir(dirOutIsoFlashes)
for i = 1:length(myFavoriteExp)
   expName = dataMatrixFlashes(myFavoriteExp(i)).expName(end-22:end);
   system(['copy ',dirIn, expName, ' ' dirOutIsoFlashes, expName])
end



%% high iso whisk
[myFavoriteExp] = findMyExpMulti(dataMatrixFlashes, [], 'iso', 1.2, whiskIndex);
dirOutIsoWhisk = [dirOut, 'Iso_whisk\'];

mkdir(dirOutIsoWhisk)
for i = 1:length(myFavoriteExp)
   expName = dataMatrixFlashes(myFavoriteExp(i)).expName(end-22:end);
   system(['copy ',dirIn, expName, ' ' dirOutIsoWhisk, expName])
end


%% high prop flashes
[myFavoriteExp] = findMyExpMulti(dataMatrixFlashes, [], 'prop', 35, flashIndex);

dirOutPropFlashes = [dirOut, 'Prop_flashes\'];

mkdir(dirOutPropFlashes)
for i = 1:length(myFavoriteExp)
   expName = dataMatrixFlashes(myFavoriteExp(i)).expName(end-22:end);
   system(['copy ',dirIn, expName, ' ' dirOutPropFlashes, expName])
end



%% high prop whisk
[myFavoriteExp] = findMyExpMulti(dataMatrixFlashes, [], 'prop', 35, whiskIndex);

dirOutPropWhisk = [dirOut, 'Prop_whisk\'];

mkdir(dirOutPropWhisk)
for i = 1:length(myFavoriteExp)
   expName = dataMatrixFlashes(myFavoriteExp(i)).expName(end-22:end);
   system(['copy ',dirIn, expName, ' ' dirOutPropWhisk, expName])
end
