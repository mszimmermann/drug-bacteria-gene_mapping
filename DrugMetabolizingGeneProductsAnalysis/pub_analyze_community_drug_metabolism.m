%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Pipeline pub_analyze_community_drug_metabolism.m
% analysis of community profiles (drugs and metabolites)
% calculate normalized intensity profiles over time and 
% drug and metabolite conversion slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables;
% takes as input infolder_DataCommunitiesDrugMetabolism = 'Data\';
% intensityNoise = 5000; % noise intensity level 
% makes an output in outfile_table_community_drug_metabolism';

fileNames = {[infolder_DataCommunitiesDrugMetabolism 'MZ012D-Drugs-'] ...
            [infolder_DataCommunitiesDrugMetabolism 'MZ012D-Metabolites-']};
        
AssayedODs = {'P1'; 'P2'}; % two plates
TimePoints = {'T00' 'T01' 'T02' 'T06' 'T08' 'T12' 'T18' 'T24'};

for fi = 1:length(fileNames)
    dataMatrix_plate = cell(size(AssayedODs));
    dataTime_plate = cell(size(AssayedODs));
    dataSamples_plate = cell(size(AssayedODs));
    for odi = 1:length(AssayedODs)
        clear RawData
        % import peak integrations from csv files, RT and AUC
        for i = 1:length(TimePoints)
               TempStruct = ReadMixedTxt([ fileNames{fi} AssayedODs{odi} '-' TimePoints{i} '_filtered.csv'], ',');

               RawData.(TimePoints{i}) = struct;
               TempStruct.IntensitiesRaw( TempStruct.IntensitiesRaw < intensityNoise)=intensityNoise; %put detection limit to 1000
               RawData.(TimePoints{i}) = TempStruct;
        end
  
        % correct by internal standard
        IS = {'IS_YOH', 'IS_CAF', 'IS_SUL'};% internal standard with which to correct intensities
        for i = 1:length(TimePoints)
            IdxIS = ismember((cellfun(@(x) x(1:6), RawData.(TimePoints{i}).Compounds, 'UniformOutput', false)), IS);
            TempMatrix = RawData.(TimePoints{i}).IntensitiesRaw;
            TempIS = TempMatrix(:,IdxIS); %intensities of internal std of this time point
            MeanIS = mean(TempIS);
            CorrectionMatrix = ones(size(TempIS));
            for p = 1:length(MeanIS)
                CorrectionMatrix(:,p)= TempIS(:,p)./MeanIS(p); %calcuate fold change from the mean
            end
            Correction = mean(CorrectionMatrix,2);
            for p = 1: length(Correction)
                TempMatrix(p,:)= TempMatrix(p,:)./Correction(p);
            end
            % eliminate samples with bad injections
            TempInjectionSum = sum(TempMatrix,2);
            TempMean = mean(TempInjectionSum);
            Outliers = find(TempInjectionSum < (TempMean/2));
            TempMatrix(Outliers, :)=NaN;

            RawData.(TimePoints{i}).IntensitiesCorrected = TempMatrix;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create matrix of community and time versus compounds
        % select compounds
        compoundNames = RawData.T01.Compounds';

        compoundsDrugs = compoundNames(cellfun(@(x) ~contains(x, 'IS_'), compoundNames));
        compoundsIDX = find(cellfun(@(x) ~contains(x, 'IS_'), compoundNames));

        dataMatrix = zeros(length(RawData.T00.SampleNames)*length(fields(RawData)),...
                           length(compoundsIDX));
        dataTime = zeros(length(RawData.T00.SampleNames)*length(fields(RawData)),1);
        dataSamples = cell(length(RawData.T00.SampleNames)*length(fields(RawData)),1);
        idx = 1;
        for j=1:length(fields(RawData))
            TempMatrix = RawData.(TimePoints{j}).IntensitiesCorrected;
            TempTime = cellfun(@(x) str2double(x(end-3:end-2)), RawData.(TimePoints{j}).SampleNames);
            TempSamples = cellfun(@(x) strsplit(x, '-'), RawData.(TimePoints{j}).SampleNames, 'unif', 0);
            TempSamples = cellfun(@(x) x{3}, TempSamples, 'unif', 0);
            TempCompounds = RawData.(TimePoints{j}).Compounds;
            % add temporary elements to the global data matrix
            [~, ~, cpdidx] = intersect(compoundsDrugs, TempCompounds, 'stable');
            dataMatrix(idx:idx+size(TempMatrix,1)-1,:) = TempMatrix(:, cpdidx);
            dataTime(idx:idx+size(TempMatrix,1)-1) = TempTime;
            dataSamples(idx:idx+size(TempMatrix,1)-1) = TempSamples;
            idx = idx+size(TempMatrix,1);
        end

        dataMatrix_plate{odi} = dataMatrix;
        dataTime_plate{odi} = dataTime;
        dataSamples_plate{odi} = dataSamples;
    end
    
    if contains(lower(fileNames{fi}), 'mz012d-drugs-')
        dataMatrixDRUG = cell2mat(dataMatrix_plate);
        dataTimeDRUG = cell2mat(dataTime_plate);
        dataSamplesDRUG = vertcat(dataSamples_plate{:});
        compoundNamesDRUG = compoundsDrugs;
    elseif contains(lower(fileNames{fi}), 'mz012d-metabolites-')
        dataMatrixMET = cell2mat(dataMatrix_plate);
        dataTimeMET = cell2mat(dataTime_plate);
        dataSamplesMET = vertcat(dataSamples_plate{:});
        compoundNamesMET = compoundsDrugs;
    end
end
clear dataMatrix dataMatrix_plate dataTime dataTime_plate 
clear dataSamples dataSamples_plate compoundNames compoundsDrugs compoundsIDX
clear CorrectionMatrix TempIS MeanIS p Correction TempMatrix
clear TempInjectionSum TempMean TempSTD Outlayers RawData TimePoints
clear TempCompounds TempSamples TempSTD TempStruct IdsIS TempTime odi fi IS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate normalized intensities meand and std and
% the slopes of drug conversion and metabolite production
for fi=1:length(fileNames)
    if fi==1
        dataMatrix = dataMatrixDRUG;
        dataTime = dataTimeDRUG;
        dataSamples = dataSamplesDRUG;
        compoundNames = compoundNamesDRUG;
        compoundIDX = 1:length(compoundNames);
        slopeDir = 1;
    elseif fi==2
        dataMatrix = dataMatrixMET;
        dataTime = dataTimeMET;
        dataSamples = dataSamplesMET;
        compoundNames = compoundNamesMET;
        compoundIDX = 1:length(compoundNames);
        slopeDir = 1;
        
    end
    dataTime_unique = unique(dataTime);
    dataCTR = cellfun(@(x) contains(x, 'CTR'), dataSamples);
    % extract unique time and samples
    dataSamplesType = cellfun(@(x) x(1:end-1), dataSamples, 'unif', 0);
    dataSamplesType_unique = unique(dataSamplesType);
    dataMatrixFC = zeros(length(dataSamplesType_unique), length(compoundNames)*length(dataTime_unique));
    dataMatrixSTD = zeros(length(dataSamplesType_unique), length(compoundNames)*length(dataTime_unique));
    dataMatrixFCcompound = cell(1, length(compoundNames)*length(dataTime_unique));
    dataMatrixFCtime = zeros(1, length(compoundNames)*length(dataTime_unique));

    dataSlopematrix = zeros(length(dataSamplesType_unique), length(compoundNames));
    fcidx = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each compound, calculate normalized intensity and conversion slope

            
    for i = 1:length(compoundIDX)  
        cmpd_i = compoundIDX(i);
        % calculate STD of data at each time point and determine the time
        % with max STD
        curSTDtime = zeros(length(dataTime_unique),1);
        for ti = 1:length(dataTime_unique)
            curSTDtime(ti) = std(dataMatrix(dataTime==dataTime_unique(ti) &...
                                ~dataCTR,...
                                 cmpd_i));
        end
        %determine the time with max STD
        [~, maxstd] = max(curSTDtime);
        
        % for each sample, calculate mean and std normalized intensity per time 
        % and conversion slope of the current compound cmpd_i
        for sample_i=1:length(dataSamplesType_unique)
            curData = dataMatrix(ismember(dataSamplesType, dataSamplesType_unique(sample_i)),...
                                 cmpd_i);
            curTime = dataTime(ismember(dataSamplesType, dataSamplesType_unique(sample_i)));
            curTime_unique = unique(curTime);
            
            curCTRL = dataMatrix(cellfun(@(x) contains(x, 'CTR'),dataSamplesType),cmpd_i);
            curCTRLtime = dataTime(cellfun(@(x) contains(x, 'CTR'),dataSamplesType));
            
            % replace low intensity with the noise level
            curData(curData<intensityNoise) = intensityNoise;
            curCTRL(curCTRL<intensityNoise) = intensityNoise;
            %calculate mean value per time point and plot as horizontal lines
            meanPerTime = zeros(length(curTime_unique),1);
            stdPerTime = zeros(length(curTime_unique),1);
            for j=1:length(curTime_unique)
                if fi==1
                    meanPerTime(j) = mean(curData(curTime==curTime_unique(j))) ./...
                                          mean(curCTRL(curCTRLtime==curTime_unique(j))) ;
                    stdPerTime(j) = meanPerTime(j)*...
                                    sqrt( (std(curData(curTime==curTime_unique(j)))/...
                                           mean(curData(curTime==curTime_unique(j))))^2 +...
                                          (std(curCTRL(curCTRLtime==curTime_unique(j)))/...
                                           mean(curCTRL(curCTRLtime==curTime_unique(j))))^2);
                elseif fi==2
                    meanPerTime(j) = mean(curData(curTime==curTime_unique(j))) ;
                    stdPerTime(j) =  std(curData(curTime==curTime_unique(j)));
                end
            end
            if fi==1
                stdPerTime = (meanPerTime/meanPerTime(1)).*...
                             sqrt( (stdPerTime./meanPerTime).^2 + ...
                                   (stdPerTime(1)/meanPerTime(1))^2);
                meanPerTime = meanPerTime/meanPerTime(1);
            end

             if fi==1
                % calculate slope to the time-point with max std across data
                [maxSlope, maxSlopeB] = calculateConversionSlope(curTime_unique,...
                                                             log2(meanPerTime),...
                                                             max(3,maxstd-1),-1);
                maxSlope = -maxSlope;
            end
            if fi==2
                % calculate slope with the sliding window of three points
                [maxSlope, maxSlopeB] = calculateConversionSlope(curTime_unique,...
                                                             meanPerTime,...
                                                             3, slopeDir);
            end
                      
            % save calculated results into summary matrices
            dataSlopematrix(sample_i, cmpd_i) = maxSlope;
            dataMatrixFC(sample_i, fcidx:fcidx+length(meanPerTime)-1) = meanPerTime;
            dataMatrixSTD(sample_i, fcidx:fcidx+length(meanPerTime)-1) = stdPerTime;            
        end
        dataMatrixFCcompound(1, fcidx:fcidx+length(meanPerTime)-1) = repmat(compoundNames(cmpd_i), size(curTime_unique));
        dataMatrixFCtime(1, fcidx:fcidx+length(meanPerTime)-1) = curTime_unique;
        fcidx = fcidx+length(meanPerTime);
    end
    if fi==1
        dataMatrixFCDRUG = dataMatrixFC;
        dataMatrixSTDDRUG = dataMatrixSTD;
        dataMatrixSlopeDRUG = dataSlopematrix;
        dataMatrixFCDRUGcompound = dataMatrixFCcompound;
        dataMatrixFCDRUGtime = dataMatrixFCtime;
    end
    if fi==2
        dataMatrixFCMET = dataMatrixFC;
        dataMatrixSTDMET = dataMatrixSTD;
        dataMatrixSlopeMET = dataSlopematrix;
        dataMatrixFCMETcompound = dataMatrixFCcompound;
        dataMatrixFCMETtime = dataMatrixFCtime;
    end
end   
%remove "Results" from the compound names
compoundNamesDRUG = cellfun(@(x) strrep(x, 'Results',''),compoundNamesDRUG,'unif',0);
compoundNamesMET = cellfun(@(x) strrep(x, 'Results',''),compoundNamesMET,'unif',0);
dataMatrixFCDRUGcompound = cellfun(@(x) strrep(x, 'Results',''),dataMatrixFCDRUGcompound,'unif',0);
dataMatrixFCMETcompound = cellfun(@(x) strrep(x, 'Results',''),dataMatrixFCMETcompound,'unif',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save normalized intensities and slopes to file
fid = fopen(outfile_table_community_drug_metabolism,'w');
fprintf(fid, 'Community');
for i=1:length(compoundNamesDRUG)
    fprintf(fid, ';Slope_%s', compoundNamesDRUG{i});
end
for i=1:length(compoundNamesMET)
    fprintf(fid, ';Slope_%s', compoundNamesMET{i});
end
for i=1:length(dataMatrixFCDRUGcompound)
    fprintf(fid, ';Mean_%s;STD_%s', dataMatrixFCDRUGcompound{i},...
                                    dataMatrixFCDRUGcompound{i});
end
for i=1:length(dataMatrixFCMETcompound)
    fprintf(fid, ';Mean_%s;STD_%s', dataMatrixFCMETcompound{i},...
                                    dataMatrixFCMETcompound{i});
end
fprintf(fid,'\n');
% printf time info to file
fprintf(fid, ';;;;;;;;');
for i=1:length(dataMatrixFCDRUGtime)
    fprintf(fid, ';%d;%d', dataMatrixFCDRUGtime(i), dataMatrixFCDRUGtime(i));
end
for i=1:length(dataMatrixFCMETtime)
    fprintf(fid, ';%d;%d', dataMatrixFCMETtime(i), dataMatrixFCMETtime(i));
end
fprintf(fid,'\n');
% print slope and intensity data to file
for i=1:length(dataSamplesType_unique)
    fprintf(fid, '%s', dataSamplesType_unique{i});
    for j=1:size(dataMatrixSlopeDRUG,2)
        fprintf(fid, ';%.3f', dataMatrixSlopeDRUG(i,j));
    end
    for j=1:size(dataMatrixSlopeMET,2)
        fprintf(fid, ';%.3f', dataMatrixSlopeMET(i,j));
    end
    for j=1:size(dataMatrixFCDRUG,2)
        fprintf(fid, ';%.3f;%.3f', dataMatrixFCDRUG(i,j),dataMatrixSTDDRUG(i,j));
    end
    for j=1:size(dataMatrixFCMET,2)
        fprintf(fid, ';%.3f;%.3f', dataMatrixFCMET(i,j),dataMatrixSTDMET(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%