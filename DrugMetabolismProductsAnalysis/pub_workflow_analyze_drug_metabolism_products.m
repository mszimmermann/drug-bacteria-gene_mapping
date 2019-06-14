%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Analyse untargeted metabolite data of the 21 drug pools 
% incubated with single bacterial species for candidate drug metabolites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_workflow_analyze_drug_metabolism_products.m
% read drug-bacteria untargeted metabolomics data
% dor each detected metabolite, calculate fold changes at t=12h and t=0h
% between specific drug-containing pools and all other pools
% save data to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % pub_workflow_analyze_drug_metabolism_products
% % takes as input infile_drug_experiment
% % takes as input infile_single_species
% infolder_drug_metabolite_data = ['DrugMetabolismProductsAnalysis' filesep 'DataFiles' filesep];
% outfile_tableS6_metaboliteIntensityFCData = ['Output' filesep 'TableS6_metabolite_intensity_and_fold_changes.csv'];
% outfile_tableS5_metabolite_filtering = ['Output' filesep 'TableS5_drug_metabolite_candidates_filtering.csv'];
% outfile_figure_pFDR_pooling_cutoff
% intensityNoise = 5000;
% fcThresholdDrug = abs(log2(0.8)); % 20% reduction
% pThreshold = 0.05;
% massThreshold = 0.002;
% RTthreshold = 0.15;
% pThresholdPooling = 10^(-6); % significance threshold for metabolite
%                              % in the drug pool vs other pools
% calculate_pThresholdPoolingFlag = 0; % whether to calculate significance
%                                      % threshold based on drug and random pooling scheme metabolite
%                                      % distributions (time consuming step)

% load species information
AssaySpecies = table2cell(readtable(infile_single_species));
% Import general experimental parameter
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);

ExperimentParameters.TimePoints = {'T00', 'T12'}; % define the timepoints to be extracted
% add control pools to the pooling scheme
addedPools = [97 98 99];
ExperimentParameters.PoolNumbers = [1:size(ExperimentParameters.PoolingScheme,2) addedPools];
ExperimentParameters.PoolingScheme = [ExperimentParameters.PoolingScheme ...
                                      zeros(size(ExperimentParameters.PoolingScheme,1), length(addedPools))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AssayConditions=cell(88,1);
idx=1;
for i = 5:88
    curAssay = ['S' repmat('0',1, 1+(i<10)), num2str(i)];
    datafile = [infolder_drug_metabolite_data curAssay '_RawDataListExportOptimized.csv'];
    if exist(datafile, 'file')==2
        AssayConditions{idx} = curAssay;
        idx=idx+1;
    else
        %sprintf('Cannot find %s', curAssay)
    end
end
% remove unised AssayConditions
AssayConditions(idx:end) = [];
clear idx curAssay

RawData = struct;

% get species names
AssayConditionsSpecies = cell(length(AssayConditions),1);
for i=1:length(AssayConditions)
    idx = cellfun(@(x) ismember(x, str2double(AssayConditions{i}(2:end))),...
                       AssaySpecies(:,1));
    AssayConditionsSpecies{i} = strjoin(AssaySpecies(idx,3:5),' ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matrices of intensity fold change between drug-containing pools
% and other pools
fcMatrix_t0_cell = cell(length(AssayConditions),1);
fdrMatrix_t0_cell= cell(length(AssayConditions),1);
fcMatrix_t12_cell = cell(length(AssayConditions),1);
fdrMatrix_t12_cell = cell(length(AssayConditions),1);
curCompounds_cell = cell(length(AssayConditions),1);
curCompoundsSpectrum_cell = cell(length(AssayConditions),1);

for d = 1:length(AssayConditions)
    datafile = [infolder_drug_metabolite_data AssayConditions{d} '_RawDataListExportOptimized.csv'];
    % Import data of one experimental condition of the single species assays (24 pools, 2 timepoints)
    curIntensitiesRaw =  Import_Intensities(datafile, 6, inf);
    curCompounds =  Import_Compounds(datafile, 6, inf);
    curSamples = Import_Injections(datafile, 5, 5)';
    curCompoundsSpectrum =  Import_CompositeSpectrum(datafile, 6, inf);
    
    %check if intensities are likely on logarithmic scale
    if nnz(curIntensitiesRaw==0)>0 && nnz(curIntensitiesRaw>100000)<1000
        curIntensitiesRaw = 10.^(curIntensitiesRaw);
    end
    % find outliers due to bad injections, if sum of all compounds is
    % lower than the mean of all 48 data sets minus three times their std
    % threshold value
    Threshold = mean(sum(curIntensitiesRaw)-(3*(std(sum(curIntensitiesRaw)))));
    curExcludedDatasets = find(sum(curIntensitiesRaw) < Threshold);
    curIntensitiesNorm = curIntensitiesRaw;
    curIntensitiesNorm(:,curExcludedDatasets) = [];
    curSamples(curExcludedDatasets) = [];
    
    % exclude ions that were unique to the excluded dataset in all
    % parameters
    curIdxExcluded = find(sum(curIntensitiesNorm,2) ==0);
    curIntensitiesNorm(curIdxExcluded,:) =  [];
    curCompounds(curIdxExcluded,:) =  [];
    curCompoundsSpectrum(curIdxExcluded,:) =  [];
   
    curIntensitiesNorm(curIntensitiesNorm==1) = NaN;
    curIntensitiesNorm = quantilenorm(curIntensitiesNorm);
    curIntensitiesNorm(isnan(curIntensitiesNorm)) = intensityNoise; %noise level
 
    totTime = cellfun(@(x) str2double(x(28:29)), curSamples);
    totDataPools = cellfun(@(x) str2double(x(strfind(x,'Pool')+(4:5))), curSamples);

    totDataMat = curIntensitiesNorm';
    % Differential analysis of ions to identify drug metabolites
    tic
    [fcMatrix_t0, fdrMatrix_t0,...
     fcMatrix_t12, fdrMatrix_t12] = workflow_DifferentialAnalysisPool(totDataMat,...
                                                                       ExperimentParameters,...
                                                                       totTime,...
                                                                       totDataPools);
    toc
    fcMatrix_t0_cell{d} = fcMatrix_t0;
    fdrMatrix_t0_cell{d} = fdrMatrix_t0;
    fcMatrix_t12_cell{d} = fcMatrix_t12;
    fdrMatrix_t12_cell{d} = fdrMatrix_t12;
    curCompounds_cell{d} = curCompounds;
    curCompoundsSpectrum_cell{d} = curCompoundsSpectrum;
    fprintf('Calculated fold changes between pools for %s (%s)\n', AssayConditions{d}, AssayConditionsSpecies{d});
    
end
clear datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded 
clear totTime totDataPools totDataMat 
clear fcMatrix_t0 fdrMatrix_t0 fcMatrix_t12 fdrMatrix_t12 fcMatrix_t0_ind fcMatrix_t12_ind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform automatic threshold calculation if
% calculate_pThresholdPoolingFlag==1
if calculate_pThresholdPoolingFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to select a p-value threshold for metabolites, 
% analyse the p-values for drugs in the pools (positive controls)
% for each drug, calculate FC and FDR
    drugFoundFC_t0 = zeros(length(ExperimentParameters.DrugMasses), length(fcMatrix_t0_cell));
    drugFoundFC_t12 = zeros(length(ExperimentParameters.DrugMasses), length(fcMatrix_t0_cell));
    drugFoundFDR_t0 = 2*ones(length(ExperimentParameters.DrugMasses), length(fcMatrix_t0_cell));
    drugFoundFDR_t12 = 2*ones(length(ExperimentParameters.DrugMasses), length(fcMatrix_t0_cell));

    for j=1:length(fcMatrix_t0_cell)
        curCompounds = curCompounds_cell{j};

        curMZ = cellfun(@(x) str2double(x(1:strfind(x, '@')-1)), curCompounds);
        curRT = cellfun(@(x) str2double(x(strfind(x, '@')+1:end)), curCompounds);

        for i=1:length(ExperimentParameters.DrugNames)

            curDrugMass = ExperimentParameters.DrugMasses(i);
            curDrugRT = ExperimentParameters.TargetedRT(i);

            drugIdx = find( (abs(curMZ-curDrugMass)<=massThreshold) &...
                            (abs(curRT-curDrugRT)<=RTthreshold) );

            if length(drugIdx)>1
                rtDiff = abs(curRT(drugIdx)-curDrugRT);
                drugIdx = drugIdx(rtDiff == min(rtDiff));
            end
            if ~isempty(drugIdx)    
                drugFoundFC_t0(i,j) = fcMatrix_t0_cell{j}(drugIdx,i);
                drugFoundFC_t12(i,j) = fcMatrix_t12_cell{j}(drugIdx,i);
                drugFoundFDR_t0(i,j) = fdrMatrix_t0_cell{j}(drugIdx,i);
                drugFoundFDR_t12(i,j) = fdrMatrix_t12_cell{j}(drugIdx,i);
            end
        end
        disp(['Done calculating pool abundances for drug ' num2str(j)])
    end
    clear rtDiff drugIdx curCompounds curMZ curRT curDrugMass curDrugRT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate pooling schemes not used originally
    pooling_patterns = [ones(1,4) zeros(1,17)];
    pooling_max = bin2dec('111100000000000000000');
    pooling_possible_pools = 1:pooling_max;
    pooling_sum = arrayfun(@(x) nnz(dec2bin(x)=='1'), pooling_possible_pools);
    pooling_possible_pools = pooling_possible_pools(pooling_sum==4);
    clear pooling_sum
    % calculate sum of the original pooling scheme and remove these from
    % possible pools
    original_pooling_scheme = binaryVectorToDecimal(ExperimentParameters.PoolingScheme(:,1:21));
    pooling_possible_pools(arrayfun(@(x) nnz(original_pooling_scheme==x)>0,...
                                          pooling_possible_pools))=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get random pooling scheme
    fcMatrix_t0_cell_random = cell(length(AssayConditions),1);
    fdrMatrix_t0_cell_random= cell(length(AssayConditions),1);
    fcMatrix_t12_cell_random = cell(length(AssayConditions),1);
    fdrMatrix_t12_cell_random = cell(length(AssayConditions),1);
    for ri = 1:randpoolnum

        random_pooling = randperm(length(pooling_possible_pools), length(ExperimentParameters.DrugNames));
        random_pooling = pooling_possible_pools(random_pooling);
        random_pooling_bin = zeros(length(random_pooling),21);
        for i=1:length(random_pooling)
            curpool = dec2bin(random_pooling(i))=='1';
            random_pooling_bin(i,22-length(curpool):end) = curpool;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ExperimentParameters.PoolingScheme = [random_pooling_bin ...
        zeros(size(ExperimentParameters.PoolingScheme,1), length(addedPools))];

        for d = 1:length(AssayConditions)
            datafile = [infolder_drug_metabolite_data AssayConditions{d} '_RawDataListExportOptimized.csv'];
            % Import data of one experimental condition of the single species assays (24 pools, 2 timepoints)
            curIntensitiesRaw =  Import_Intensities(datafile, 6, inf);
            curSamples = Import_Injections(datafile, 5, 5)';

            % find outliers due to bad injections, if sum of all compounds is
            % lower than the mean of all 48 data sets minus three times their std
            % threshold value

            Threshold = mean(sum(curIntensitiesRaw)-(3*(std(sum(curIntensitiesRaw)))));
            curExcludedDatasets = find(sum(curIntensitiesRaw) < Threshold);
            curIntensitiesNorm = curIntensitiesRaw;
            curIntensitiesNorm(:,curExcludedDatasets) = [];
            curSamples(curExcludedDatasets) = [];

            % exclude ions that were unique to the excluded dataset in all
            % parameters
            curIdxExcluded = sum(curIntensitiesNorm,2) ==0;
            curIntensitiesNorm(curIdxExcluded,:) =  [];

            curIntensitiesNorm(curIntensitiesNorm==1) = NaN;
            curIntensitiesNorm = quantilenorm(curIntensitiesNorm);
            curIntensitiesNorm(isnan(curIntensitiesNorm)) = intensityNoise; %noise level

            totTime = cellfun(@(x) str2double(x(28:29)), curSamples);
            totDataPools = cellfun(@(x) str2double(x(strfind(x,'Pool')+(4:5))), curSamples);

            totDataMat = curIntensitiesNorm';
            % Differential analysis of ions to identify drug metabolites
            tic
            [fcMatrix_t0, fdrMatrix_t0,...
             fcMatrix_t12, fdrMatrix_t12] = workflow_DifferentialAnalysisPool(totDataMat,...
                                                                               ExperimentParameters,...
                                                                               totTime,... 
                                                                               totDataPools);

            toc
            fcMatrix_t0_cell_random{d} = fcMatrix_t0;
            fdrMatrix_t0_cell_random{d} = fdrMatrix_t0;
            fcMatrix_t12_cell_random{d} = fcMatrix_t12;
            fdrMatrix_t12_cell_random{d} = fdrMatrix_t12;
            fprintf('Calculated fold changes between random pools for %s (%s)\n', AssayConditions{d}, AssayConditionsSpecies{d});
        end
    end
    clear datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
    clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded 
    clear totTime totDataPools totDataMat 
    clear fcMatrix_t0 fdrMatrix_t0 fcMatrix_t12 fdrMatrix_t12 fcMatrix_t0_ind fcMatrix_t12_ind

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make a histogram of all p-values at t12, 
    % p-values for drugs at t0 (positive controls)
    % histogram of all p-values at t12 of random pool (negative control)
    randpnum = 5000;
    allFDRt12 = zeros(randpnum*length(fdrMatrix_t12_cell)+1,size(fdrMatrix_t12_cell{1},2));
    allFDRt12rand = zeros(randpnum*length(fdrMatrix_t12_cell)+1,size(fdrMatrix_t12_cell{1},2));

    idx = 1;
    for i=1:length(fdrMatrix_t12_cell)
        randmets = randi(size(fdrMatrix_t12_cell{i},1), 1, randpnum);
        allFDRt12(idx:idx+randpnum-1,:) = fdrMatrix_t12_cell{i}(randmets,:);
        randmets = randi(size(fdrMatrix_t12_cell_random{i},1), 1, randpnum);
        allFDRt12rand(idx:idx+randpnum-1,:) = fdrMatrix_t12_cell_random{i}(randmets,:);
        idx = idx+randpnum;
    end
    allFDRt12(idx:end,:) = [];
    allFDRt12rand(idx:end,:) = [];
    % make a histogram of allFDRt12
    [N12,edges12] = histcounts(-log10(allFDRt12(1:10:end)),1000);
    edges12 = edges12(1:end-1) + (edges12(2:end)-edges12(1:end-1))/2;
    % make a histogram of real drug p-values in drug pools
    [NDrug,edgesDrug] = histcounts(-log10(drugFoundFDR_t0(:)),1000);
    edgesDrug = edgesDrug(1:end-1) + (edgesDrug(2:end)-edgesDrug(1:end-1))/2;
    [N12rand,edges12rand] = histcounts(-log10(allFDRt12rand(1:10:end)),1000);
    edges12rand = edges12rand(1:end-1) + (edges12rand(2:end)-edges12rand(1:end-1))/2;

    clear allFDRt12 allFDRt12rand
  
    % plot threeo histograms together
    figure
    plot(edges12, N12, 'k', 'LineWidth',2)
    set(gca, 'YScale', 'log')
    hold on
    plot(edgesDrug, NDrug, 'r')
    set(gca, 'YScale', 'log')
    hold on
    plot(edges12rand, N12rand, 'b', 'LineWidth',1)
    set(gca, 'YScale', 'log')

    % smooth the curve and find the beginning of the significant peak
    smoothFDRDrug = smooth(NDrug(1:20:end));
    smoothedgesDrug = edgesDrug(1:20:end);
    %plot(smoothedgesDrug, smoothFDRDrug, 'r', 'LineWidth', 2)
    [testpeaks, testloc] = findpeaks(-smoothFDRDrug(1:length(smoothFDRDrug)/2));
    hold on
    for i=1:length(testloc)
        plot([floor(smoothedgesDrug(testloc(i))), floor(smoothedgesDrug(testloc(i)))], [1, 10^10], '--')
    end
    smoothedgesDrug(testloc);
    xlim([-5 40])
    ylim([1 10^7])
    %legend({'FDR distribution, t12', 'Drug FDR distribution, t12', 'Smoothed drug FDR', 'Automatic cutoff'})
    legend({'FDR distribution, t12 (subset)',...
            'Drug FDR distribution, t0 (positive control)',...
            'Random pooling FDR (subset) (negative control)',...
            'Automatic cutoff'})
    xlabel('-log10(FDR) in drug pool vs other pools')
    ylabel('Number of metabolites')
    title('Significance of pool-specific metabolites, t12')
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
          outfile_figure_pFDR_pooling_cutoff)
    pThresholdPooling = 10^(-floor(smoothedgesDrug(testloc)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create 3D matrix with metabolites,drugs and species
% first make a list of changing metabolites
changingMets = cell(1000000,1);
changingMetsSpectrum = cell(1000000,1);
changingMetsPool = zeros(1000000,1);

idx = 1;
for i=1:length(curCompounds_cell)
    metFC = log2(fcMatrix_t0_cell{i});
    metPadj = fdrMatrix_t0_cell{i};
    curChange = sum((metFC>=fcThresholdPool) & (metPadj <= pThresholdPooling),2);
    curMets = curCompounds_cell{i}(curChange>0);
    curMetsSpectrum = curCompoundsSpectrum_cell{i}(curChange>0);
    %add only new metabolites
    [curMets, idxmets] = setdiff(curMets, changingMets(1:idx-1));
    curMetsSpectrum = curMetsSpectrum(idxmets);
    changingMets(idx:idx+length(curMets)-1) = curMets;
    changingMetsSpectrum(idx:idx+length(curMets)-1) = curMetsSpectrum;
    changingMetsPool(idx:idx+length(curMets)-1)=i;
    idx = idx+length(curMets);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add T12 metabolites
    metFC = log2(fcMatrix_t12_cell{i});
    metPadj = fdrMatrix_t12_cell{i};
    curChange = sum(metFC>=fcThresholdPool & metPadj <= pThresholdPooling,2);
    curMets = curCompounds_cell{i}(curChange>0);
    curMetsSpectrum = curCompoundsSpectrum_cell{i}(curChange>0);
    %add only new metabolites
    [curMets,idxmets] = setdiff(curMets, changingMets(1:idx-1));
    curMetsSpectrum = curMetsSpectrum(idxmets);
    changingMets(idx:idx+length(curMets)-1) = curMets;
    changingMetsSpectrum(idx:idx+length(curMets)-1) = curMetsSpectrum;
    changingMetsPool(idx:idx+length(curMets)-1)=i;
    idx = idx+length(curMets);
    fprintf('Identified metabolites in drug pools for %s (%s)\n', AssayConditions{i}, AssayConditionsSpecies{i});
end
changingMets(idx:end) = [];
changingMetsSpectrum(idx:end) = [];
changingMetsPool(idx:end) = [];
clear curMets curChange metFC metPadj
[changingMets,idx] = unique(changingMets);
changingMetsSpectrum = changingMetsSpectrum(idx);
changingMetsPool = changingMetsPool(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge metabolites that are close in mass and RT
[changingMets_merged,...
          changingMets_merged_idx,...
          changingMets_merged_idx_unique,...
          changingMets_merged_spectrum,...
          changingMets_merged_mass,...
          changingMets_merged_RT,...
          changingMets_merged_mass_delta,...
          changingMets_merged_RT_delta,...
          changingMets_merged_number] = merge_changing_metabolites(changingMets,...
                                                                   changingMetsSpectrum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get idx conversion from each ion file to the merged
curCompoundsMergedIDXconversion = zeros(length(changingMets_merged), length(curCompounds_cell));
for d = 1:length(curCompounds_cell)
    curMZ = cellfun(@(x) str2double(x(1:strfind(x, '@')-1)), curCompounds_cell{d});
    curRT = cellfun(@(x) str2double(x(strfind(x, '@')+1:end)), curCompounds_cell{d});

    tic
    for i=1:size(changingMets_merged_mass,1)

        curMetMass = changingMets_merged_mass(i);
        curMetRT = changingMets_merged_RT(i);

        metIdx = find( (abs(curMZ-curMetMass)<=max(massThreshold,...
                                                   changingMets_merged_mass_delta(i))) &...
                        (abs(curRT-curMetRT)<=max(RTthreshold,...
                                                  changingMets_merged_RT_delta(i))));

        if length(metIdx)>1
            rtDiff = abs(curRT(metIdx)-curMetRT);
            metIdx = metIdx(rtDiff == min(rtDiff));
            if length(metIdx)>1
                mzDiff = abs(curMZ(metIdx)-curMetMass);
                metIdx = metIdx(mzDiff == min(mzDiff));
            end
        end
        if ~isempty(metIdx)    
            curCompoundsMergedIDXconversion(i,d) = metIdx;
        end
    end
    toc
    fprintf('Merged metabolites for %s (%s)\n', AssayConditions{d}, AssayConditionsSpecies{d});
end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create metabolite-drug matrix to save in how many species metabolite is
% changing in each drug pool
changingMets_drug_matrix = zeros(length(changingMets_merged), length(ExperimentParameters.DrugNames));
for i=1:length(curCompounds_cell)
    curIDXmerged = curCompoundsMergedIDXconversion(:,i)~=0;
    curIDXfull = curCompoundsMergedIDXconversion(curIDXmerged,i);
  
    % check whether metabolites occur at t=0 in the drug pools
    metFC = log2(fcMatrix_t0_cell{i}(curIDXfull,:));
    metPadj = fdrMatrix_t0_cell{i}(curIDXfull,:);
    % add flags to changingMets_drug_matrix if metabolites are changing
    % in this species among drugs
    changingMets_drug_matrix(curIDXmerged,:) = ...
        changingMets_drug_matrix(curIDXmerged,:) + ...
        ((metFC>=fcThresholdPool) & (metPadj <= pThresholdPooling));
    
    % check whethe metabolites occur at t=12 in the drug pools
    metFC = log2(fcMatrix_t12_cell{i}(curIDXfull,:));
    metPadj = fdrMatrix_t12_cell{i}(curIDXfull,:);
    % add flags to changingMets_drug_matrix if metabolites are changing
    % in this species among drugs
    changingMets_drug_matrix(curIDXmerged,:) = ...
        changingMets_drug_matrix(curIDXmerged,:) + ...
        ((metFC>=fcThresholdPool) & (metPadj <= pThresholdPooling));
    fprintf('Calculated metabolite-drug matrix for %s (%s)\n', AssayConditions{i}, AssayConditionsSpecies{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a matrix of metabolite-drug indeces
drug_mets = zeros(nnz(changingMets_drug_matrix),2);
idx = 1;
for i=1:size(changingMets_drug_matrix,1)
    curidx = find(changingMets_drug_matrix(i,:));
    drug_mets(idx:idx+length(curidx)-1,:) = [i*ones(size(curidx')) curidx'];
    idx = idx+length(curidx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for all drug metabolites, create tables of
% intensities at t=0h and t=12h, 
% fold changes in the drug vs non-drug pools at t=0 and t=12, 
% and fold changes between t=12h and t=0h
metFoundFC_t0 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFC_t12 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFDR_t0 = 2*ones(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFDR_t12 = 2*ones(size(drug_mets,1), length(fcMatrix_t0_cell));
metIntensity_t0 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
metIntensity_t12 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
stdIntensity_t0 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
stdIntensity_t12 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFC_t12t0 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFDR_t12t0 = 2*ones(size(drug_mets,1), length(fcMatrix_t0_cell));
metFoundFCSTD_t12t0 = zeros(size(drug_mets,1), length(fcMatrix_t0_cell));

% save individual intensities for dexamethasone metabolite for plotting
selectedMetIDX = find(abs(changingMets_merged_mass - 332.179) < massThreshold);
met_DexInt_t0 = zeros(length(AssayConditions), 4);
met_DexInt_t12 = zeros(length(AssayConditions), 4);


for d = 1:length(AssayConditions)
    datafile = [infolder_drug_metabolite_data AssayConditions{d} '_RawDataListExportOptimized.csv'];
    % Import data of one experimental condition of the single species assays (24 pools, 2 timepoints)
    curIntensitiesRaw =  Import_Intensities(datafile, 6, inf);
    curSamples = Import_Injections(datafile, 5, 5)';
    
    %check if intensities are likely on logarithmic scale
    if nnz(curIntensitiesRaw==0)>0 && nnz(curIntensitiesRaw>100000)<1000
        curIntensitiesRaw = 10.^(curIntensitiesRaw);
    end
    % find outliers due to bad injections, if sum of all compounds is
    % lower than the mean of all 48 data sets minus three times their std
    % threshold value
    
    Threshold = mean(sum(curIntensitiesRaw)-(3*(std(sum(curIntensitiesRaw)))));
    curExcludedDatasets = find(sum(curIntensitiesRaw) < Threshold);
    curIntensitiesNorm = curIntensitiesRaw;
    curIntensitiesNorm(:,curExcludedDatasets) = [];
    curSamples(curExcludedDatasets) = [];
    
    % exclude ions that were unique to the excluded dataset in all
    % parameters
    curIdxExcluded = sum(curIntensitiesNorm,2) == 0;
    curIntensitiesNorm(curIdxExcluded,:) =  [];
   
    curIntensitiesNorm(curIntensitiesNorm==1) = NaN;
    curIntensitiesNorm = quantilenorm(curIntensitiesNorm);
    curIntensitiesNorm(isnan(curIntensitiesNorm)) = intensityNoise; %noise level
 
    totTime = cellfun(@(x) str2double(x(28:29)), curSamples);
    totDataPools = cellfun(@(x) str2double(x(strfind(x,'Pool')+(4:5))), curSamples);

    totDataMat = curIntensitiesNorm';
    % Differential analysis of ions to identify drug metabolites
    tic
    [fcMatrix_t0t12, fdrMatrix_t0t12, stdMatrix_t0t12,...
     intMatrix_t0mean, intMatrix_t0std,...
     intMatrix_t12mean, intMatrix_t12std,...
     rawIntensity_t0,rawIntensity_t12] = workflow_DifferentialAnalysisT0T12(totDataMat,...
                                                                           ExperimentParameters,...
                                                                           totTime,...
                                                                           totDataPools);                                                                      
    toc

    % select metabolites that are changing
    curIDXmerged = find(curCompoundsMergedIDXconversion(:,d)~=0);
    curIDXfull = curCompoundsMergedIDXconversion(curIDXmerged,d);
    
    fcMatrix_t0t12 = fcMatrix_t0t12(curIDXfull,:);
    fdrMatrix_t0t12 = fdrMatrix_t0t12(curIDXfull,:);
    stdMatrix_t0t12 = stdMatrix_t0t12(curIDXfull,:);
    intMatrix_t0mean = intMatrix_t0mean(curIDXfull,:);
    intMatrix_t0std = intMatrix_t0std(curIDXfull,:);
    intMatrix_t12mean = intMatrix_t12mean(curIDXfull,:);
    intMatrix_t12std = intMatrix_t12std(curIDXfull,:);
    fcMatrix_t0 = fcMatrix_t0_cell{d}(curIDXfull,:);
    fdrMatrix_t0 = fdrMatrix_t0_cell{d}(curIDXfull,:);
    fcMatrix_t12 = fcMatrix_t12_cell{d}(curIDXfull,:);
    fdrMatrix_t12 = fdrMatrix_t12_cell{d}(curIDXfull,:);
    
    for j=1:length(rawIntensity_t0)
        rawIntensity_t0{j} = rawIntensity_t0{j}(:,curIDXfull);
        rawIntensity_t12{j} = rawIntensity_t12{j}(:,curIDXfull);
    end
            
    for j=1:length(curIDXmerged)
        curidx = find(drug_mets(:,1)==curIDXmerged(j));
        curdrugidx = drug_mets(curidx,2);
        for k=1:length(curdrugidx)
            metFoundFC_t0(curidx(k),d) = fcMatrix_t0(j,curdrugidx(k));
            metFoundFC_t12(curidx(k),d) = fcMatrix_t12(j,curdrugidx(k));
            metFoundFDR_t0(curidx(k),d) = fdrMatrix_t0(j,curdrugidx(k));
            metFoundFDR_t12(curidx(k),d) = fdrMatrix_t12(j,curdrugidx(k));
            metIntensity_t0(curidx(k),d) = intMatrix_t0mean(j,curdrugidx(k));
            metIntensity_t12(curidx(k),d) = intMatrix_t12mean(j,curdrugidx(k));
            stdIntensity_t0(curidx(k),d) = intMatrix_t0std(j,curdrugidx(k));
            stdIntensity_t12(curidx(k),d) = intMatrix_t12std(j,curdrugidx(k));
            metFoundFC_t12t0(curidx(k),d) = fcMatrix_t0t12(j,curdrugidx(k));
            metFoundFDR_t12t0(curidx(k),d) = fdrMatrix_t0t12(j,curdrugidx(k));
            metFoundFCSTD_t12t0(curidx(k),d) = stdMatrix_t0t12(j,curdrugidx(k));
            % only save individual intensities for C.scindens for plotting
            if drug_mets(curidx(k),1) == selectedMetIDX
                met_DexInt_t0(d,:) = rawIntensity_t0{curdrugidx(k)}(:,j);
                met_DexInt_t12(d,:) = rawIntensity_t12{curdrugidx(k)}(:,j);
            end
        end
    end
    fprintf('Calculated selected metabolite fold changes between t=12h and t=0h for %s (%s)\n', AssayConditions{d}, AssayConditionsSpecies{d});

end
clear  datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded 
clear totTime totDataPools totDataMat  totDataRT totDataSampleNames
clear totCompounds totCompounds_abbr totCtrl
clear TempResultsMetabolites TempStruct
clear fcMatrix_t0t12 fdrMatrix_t0t12 stdMatrix_t0t12
clear intMatrix_t0mean intMatrix_t0std intMatrix_t12mean intMatrix_t12std
clear fcMatrix_t0 fdrMatrix_t0 fcMatrix_t12 fdrMatrix_t12
clear rawIntensity_t0 rawIntensity_t12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save individual dexamethasone metabolite info to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the file for writing
fid = fopen(outfile_table_dexamethasone_metabolite_intensityT12, 'w');
% print headers
fprintf(fid, 'ParentDrug;MZ;RT;MZdelta;Species;A;B;C;D\n');
for i=1:size(met_DexInt_t12,1)
    fprintf(fid,'Dexamethasone;%.3f;%.3f;%.3f;%s;%.3f;%.3f;%.3f;%.3f\n',...
        changingMets_merged_mass(selectedMetIDX),...
        changingMets_merged_RT(selectedMetIDX),...
        changingMets_merged_mass_delta(selectedMetIDX),...
        AssayConditionsSpecies{i},...
        met_DexInt_t12(i,:));
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save metabolite info to file
% metabolite intensities and fold changes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the file for writing
fid = fopen(outfile_tableS6_metaboliteIntensityFCData, 'w');
% print headers
fprintf(fid, 'Index;ParentDrug;MZ;RT;MZdelta');
for i=1:length(AssayConditionsSpecies)
    fprintf(fid,';%s;;;;;;;;;;',...
       AssayConditionsSpecies{i});
end
fprintf(fid, '\nIndex;ParentDrug;MZ;RT;MZdelta');
for i=1:length(AssayConditionsSpecies)
    fprintf(fid,';FC (log2) drug vs no-drug pool t=0h %s;FDR drug vs no-drug pool t=0h %s;FC (log2) drug vs no-drug pool t=12h %s;FDR drug vs no-drug pool t=12h %s;Intensity mean t=0h %s;Intensity STD t=0h %s;Intensity mean t=12h %s;Intensity STD t=12h %s;FC t=12h vs t=0h (log2) %s;FC STD t=12h vs t=0h %s;p(FDR) t=12h vs t=0h %s',...
        AssayConditionsSpecies{i},AssayConditionsSpecies{i},...
        AssayConditionsSpecies{i},AssayConditionsSpecies{i},...
        AssayConditionsSpecies{i},AssayConditionsSpecies{i},...
        AssayConditionsSpecies{i},AssayConditionsSpecies{i},...
        AssayConditionsSpecies{i},AssayConditionsSpecies{i},AssayConditionsSpecies{i});
end
fprintf(fid, '\n');
for i=1:size(drug_mets,1)
    fprintf(fid,'%d;%s;%.3f;%.3f;%.3f',...
        i,...
        ExperimentParameters.DrugNames{drug_mets(i,2)},...
        changingMets_merged_mass(drug_mets(i,1)),...
        changingMets_merged_RT(drug_mets(i,1)),...
        changingMets_merged_mass(drug_mets(i,1)) - ExperimentParameters.DrugMasses(drug_mets(i,2)));
    for j=1:length(AssayConditionsSpecies)
        log2t0 = log2(metFoundFC_t0(i,j));
        log2t12 = log2(metFoundFC_t12(i,j));
        log2t0t12 = log2(metFoundFC_t12t0(i,j));
        log2t0(isinf(log2t0)) = nan;
        log2t12(isinf(log2t12)) = nan;
        log2t0t12(isinf(log2t0t12)) = nan;
        fprintf(fid,';%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f',...
                log2t0,...
                metFoundFDR_t0(i,j),...
                log2t12,...
                metFoundFDR_t12(i,j),...
                metIntensity_t0(i,j),...
                stdIntensity_t0(i,j),...
                metIntensity_t12(i,j),...
                stdIntensity_t12(i,j),...
                log2t0t12 ,...
                metFoundFCSTD_t12t0(i,j)/(metFoundFC_t12t0(i,j)*log(2)),...
                metFoundFDR_t12t0(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate occurance of each metabolite in the list
drugmet_occurence = zeros(size(drug_mets,1), 1);
for i=1:length(drug_mets)
    drugmet_occurence(i) = nnz(drug_mets(:,1) == drug_mets(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read clustergram info to make a flag of whether the drug was metabolized
% by at least one bacterium
drugAnalysisClustergram = readtable(outfile_table_FCclustergram,...
                                   'HeaderLines',1); %skip first line);
% get the names of consumed drugs from the clustergram
clustDrugNames = drugAnalysisClustergram.DrugName;
clustDrugNames = cellfun(@(x) strtrim(x), clustDrugNames, 'unif',0);
ExperimentParameters.DrugNames = cellfun(@(x) strtrim(x), ExperimentParameters.DrugNames, 'unif',0);
% make drug metabolism flag
drug_metabolized = ismember(ExperimentParameters.DrugNames(drug_mets(:,2)),...
                            clustDrugNames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate in how many species each metabolite was changing between t=12 and t=0
drugmet_change12_0 = sum( (log2(metFoundFC_t12t0) >= fcThresholdMetT12T0Pool) &...
                          (metFoundFDR_t12t0 <= pThreshold), 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mass delta between metabolite and drugs
drugmets_drugmass_delta = changingMets_merged_mass(drug_mets(:,1)) - ...
                          ExperimentParameters.DrugMasses(drug_mets(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLAG metabolites with the same mass as parent drugs
drugmet_parentDrugs = abs(changingMets_merged_mass(drug_mets(:,1)) - ...
                          ExperimentParameters.DrugMasses(drug_mets(:,2)))<=massThreshold &...
                      abs(changingMets_merged_RT(drug_mets(:,1)) - ...
                          ExperimentParameters.TargetedRT(drug_mets(:,2))) <= RTthreshold;
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% join mass deltas of good mets whithin +/-smoothMZ interval
%smooth drugmets_drugmass_delta_trunc
drugmets_drugmass_delta_trunc = drugmets_drugmass_delta;

% smoothing window is defined in Drug_bacteria_gene_mapping_variables
% smoothMZdeltas = 0.002;
for i=1:length(drugmets_drugmass_delta_trunc)
    similaridx = find(abs(drugmets_drugmass_delta_trunc-drugmets_drugmass_delta_trunc(i))<smoothMZdeltas);
    maxidx = max(drugmets_drugmass_delta_trunc(similaridx));
    similaridx = unique([similaridx; find(abs(drugmets_drugmass_delta_trunc-maxidx)<smoothMZdeltas)]);
    drugmets_drugmass_delta_trunc(similaridx) = mean(drugmets_drugmass_delta_trunc(similaridx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(outfile_tableS5_metabolite_filtering, 'w');
fprintf(fid, ['Index;MZ;RT;ParentDrug;GoodFilter;FilterNoisyT0;FilterNoisyT12;'...
              'FilterDifferentMassDefect;FilterLikelyDrugFragment;'...
              'FilterLikelyDrugMetaboliteFragment;PotentialFragments;'...
              'IonSpectrum;LeadingMZ;DrugMassFlag;DrugMassDelta;DrugMassDeltaSmoothed;'...
              'DrugConsumedFlag;NumberOfDrugs;NumberOfIncreasedT12vsT0\n']);

for drugidx = 1:length(ExperimentParameters.DrugNames)
    curmets = drug_mets(:,2)==drugidx;
    % create current drug metabolite infor for filtering
    curdrug_mets = drug_mets(curmets,1);     
    curdrug_mets(:,2:3) = [changingMets_merged_mass(curdrug_mets(:,1)) changingMets_merged_RT(curdrug_mets(:,1))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter candidates based on pooling comparisons
    [metFilter, metPotentialFragments] = filterDrugMetaboliteHits(metIntensity_t0(curmets,:),...
                                      metIntensity_t12(curmets,:),...
                                      metFoundFC_t0(curmets,:),...
                                      metFoundFC_t12(curmets,:),...
                                      curdrug_mets,...
                                      ExperimentParameters.DrugMasses(drugidx),...
                                      ExperimentParameters.TargetedRT(drugidx));
    % leave only good candidates and make a clustergram
    good_candidates = sum(metFilter,2)==0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write to file all candidate metabolites and flags
    curmets = find(curmets); 
    for i=1:size(curdrug_mets,1)
        % get metabolite specra and leading ion
        curspectra = changingMets_merged_spectrum(curdrug_mets(i,1));
        maxMZ = parseMaxCompositeSpectra(curspectra);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid,'%d;%.3f;%.3f;%s;%d;%d;%d;%d;%d;%d',...
                    curmets(i),...
                    curdrug_mets(i,2:3),...
                    ExperimentParameters.DrugNames{drugidx},...
                    good_candidates(i),...
                    metFilter(i,:));
        fprintf(fid,';');
        for j=1:nnz(metPotentialFragments(i,:))
            fprintf(fid,'[[%d],%.3f,%.3f]',curmets(metPotentialFragments(i,j)),...
                                           curdrug_mets(metPotentialFragments(i,j),2:3));
        end
        fprintf(fid,';%s;%.3f;%d;%.3f;%.3f;%d;%d;%d',...
                    curspectra{1},...
                    maxMZ,...
                    drugmet_parentDrugs(curmets(i)),...
                    drugmets_drugmass_delta(curmets(i)),...
                    drugmets_drugmass_delta_trunc(curmets(i)),...
                    drug_metabolized(curmets(i)),...
                    drugmet_occurence(curmets(i)),...
                    drugmet_change12_0(curmets(i)));
        fprintf(fid, '\n');
    end
    disp(['Done filtering metabolites for drug' num2str(drugidx)])
end
fclose(fid);              

