%do some FSCF stats
load('Z:\code\Jennifer_code\SpikeSorting\dataFiles\poissonSurrogateData\poissonSurrogateDataAwakeFSCF.mat')


%% 
numCells = cellfun(@numel, realMIs);
numCells = sum(numCells(1,:));
kruskalInput = zeros(numCells, 6);
%for each condition
for cond = 1:6
    %take all MIs
    kruskalInput(:, cond) = cell2mat(realMIs(cond, :))';
end

%kruskal wallis the thing
 [P,ANOVATAB] = kruskalwallis(kruskalInput)
%add labels
xticklabels({"2%", "11%", "25%", "44%", "75%", "100%"});
xlabel('Contrast Level')
ylabel('MI')

%% 
numCells = cellfun(@numel, realMIs);
numCells = sum(numCells(1,:));
allSigs = zeros(numCells, 6);
%for each condition
for cond = 1:6
    %take all MIs
    allSigs(:, cond) = cell2mat(significance(cond, :))';
end

%% check out significance by bootstrap?
low95 = zeros(1,6);
high95 = zeros(1,6);
allMeans = zeros(1,6);
bootDists = zeros(100, 6);
%for each significance value
for cond = 1:6
    [bootRes, bootDists(:, cond)] = bootci(100, @nanmean, allSigs(:, cond));
    low95(cond) = bootRes(1);
    high95(cond) = bootRes(2);
    allMeans(cond) = nanmean(allSigs(:, cond));
end

figure('color', 'white', 'renderer', 'painters')
cons = [2, 11, 25, 44, 75, 100];
errorbar(cons, allMeans, allMeans - low95, high95 - allMeans, 'o-')
xlabel('Contrast')
ylabel('Fraction Entrained')
%% to excel
figureSup15 = [allMeans; low95; high95];

figureSup15 = array2table(figureSup15);

filename = 'Z:\code\Jennifer_code\SpikeSorting\AwakeExcel\figureSup15FSCF.xlsx';
writetable(figureSup15,filename,'Sheet',1,'Range','A1')
%% 
[P,ANOVATAB] = kruskalwallis(bootDists)

%% t test the boot dists?
%for each contrast
tVals = nan(6);
p_val = nan(6);
stds = nan(6, 1);
for cond1 = 1:6
    stds(cond1) = nanstd(bootDists(:, cond1));
end
for cond1 = 1:6
    %for every other contrast
    for cond2 = cond1+1:6
        mean1 = nanmean(allSigs(:, cond1));
        dist2 = bootDists(:, cond2); 
        
        
        %tVals(cond1, cond2) = (nanmean(dist1) - nanmean(dist2))/sqrt((nanstd(dist1)^2)/100 + (nanstd(dist2)^2)/100);
        %[~, ~, ~, tValsMatlab(cond1, cond2)] = ttest2(dist1', dist2');
        %p_val(cond1,cond2) = normcdf(mean1, nanmean(dist2), nanstd(dist2));
        tVals(cond1, cond2) = (allMeans(cond2) - allMeans(cond1))/stds(cond2);
    end
end

disp(tVals);
pVals = (1 - tcdf(abs(tVals), 99))*15*2;

for i = 6:-1:1
    for j = i-1:-1:1
        disp(pVals(j, i))
    end
end
%% quick plotting


for cond = 1:6
    scatter(cons(cond) * ones(503, 1), kruskalInput(:, cond));
    hold on
end

xlabel('contrast')
ylabel('MI')

%% 
% figure
% violinplot(kruskalInput);
% 
% xlabel('contrast')
% ylabel('MI')
% %% ranksum
% %for each contrast
% p = nan(6);
% for cond1 = 1:6
%     %for every other contrast
%     for cond2 = cond1+1:6
%         p(cond1, cond2) = ranksum(kruskalInput(:, cond1), kruskalInput(:, cond2));
%     end
% end
%% cross checck how many good fscf cells there are
dataPath = 'D:\adeeti\KiloSort_results\FSCFAwakeOnly\'; %CHANGE THIS 
experimentFolders = dir([dataPath, 'Output*']);
totCells = 0;
numSigCell = zeros(6, 11);
numCell = zeros(6, 11);
%for each experiment
for expInd = 1:11
    %get the location of each cell
    [~, ~, unitLayerCat, ~, ~] = getUsefulVars(experimentFolders(expInd).name(8:end), 'Awake', 'FSCF');
    inCortex = find(~isnan(unitLayerCat));
    totCells = totCells + numel(inCortex);
    numCell(:, expInd) =  numel(inCortex);
    %cell indices that are in the cortex only pls
    %for each intensity
    for conInd = 1:6
        %count significant ones in cortex only
        numSigCell(conInd, expInd) = nnz(significance{conInd, expInd}(inCortex) == 1);
    end
end
totCells
%% tmcomptest
tmInput(1:6, 1) = numSignificant; %reshaping for function
tmInput(1:6, 2) = sum(numLayer(1,:)).*ones(6,1); %reshaping for function
% load('poissonSurrogateDataAwakeLED.mat', 'numSignificant')
% %tmInput(7, :) = sum(numSignificant, 2)'; %add LED
% %tmInput(7, 2) = tmInput(7, 1)  + tmInput(7, 2); %make it a ratio

tmcomptest(tmInput, .05)

%% 
thing = numSigCell./numCell;

figure
for mouse = 1:11
    scatter([2 11 25 44 75 100], thing(:, mouse));
    hold on
end

xlabel('contrast')
ylabel('proportion entrained')