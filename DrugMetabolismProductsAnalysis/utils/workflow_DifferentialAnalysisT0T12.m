function [fcMatrix_t0t12, fdrMatrix_t0t12, stdMatrix_t0t12,...
          rawIntensityMean_t0,rawIntensitySTD_t0,...
          rawIntensityMean_t12,rawIntensitySTD_t12,...
          rawIntensity_t0,rawIntensity_t12] = workflow_DifferentialAnalysisT0T12(totDataMat,...
                                                                           ExperimentParameters,...
                                                                           totTime,...
                                                                           totDataPools)
% perform differential analysis of metabolite data in drug pool vs all
% other pools at t=0 and t=12

fcMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
pMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
stdMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));

rawIntensityMean_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensityMean_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensitySTD_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensitySTD_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));

rawIntensity_t0 = cell(1,length(ExperimentParameters.DrugNames));
rawIntensity_t12 = cell(1,length(ExperimentParameters.DrugNames));

for j=1:length(ExperimentParameters.DrugNames)
    IdxPool = ExperimentParameters.PoolNumbers(ExperimentParameters.PoolingScheme(j,:)==1);
   
    curData_t0 = totDataMat(totTime==0 & ismember(totDataPools, IdxPool),:);
    curData_t12 = totDataMat(totTime==12 & ismember(totDataPools, IdxPool),:);
    
    rawIntensityMean_t0(:,j) = nanmean(curData_t0);
    rawIntensitySTD_t0(:,j) = nanstd(curData_t0);
    rawIntensityMean_t12(:,j) = nanmean(curData_t12);
    rawIntensitySTD_t12(:,j) = nanstd(curData_t12);
    
    fcMatrix_t0t12(:,j) = (nanmean(curData_t12) ./...
                           nanmean(curData_t0))';
    [~, pMatrix_t0t12(:,j)] = ttest2(curData_t12,...
                                  curData_t0,...
                                  'VarType', 'equal');
    stdMatrix_t0t12(:,j) = abs(fcMatrix_t0t12(:,j)).*...
                               sqrt( (nanstd(curData_t12)./nanmean(curData_t12)).^2 + ...
                                     (nanstd(curData_t0)./nanmean(curData_t0)).^2 )';
    
    rawIntensity_t0{j} = curData_t0;
    rawIntensity_t12{j} = curData_t12;
end
% adjust FDR
fdrMatrix_t0t12 = mafdr(pMatrix_t0t12(:), 'BHFDR', 1);
fdrMatrix_t0t12 = reshape(fdrMatrix_t0t12, size(fcMatrix_t0t12));