
function [fcMatrix_t0, fdrMatrix_t0,...
          fcMatrix_t12, fdrMatrix_t12] = workflow_DifferentialAnalysisPool(totDataMat,...
                                                                           ExperimentParameters,...
                                                                           totTime,...
                                                                           totDataPools)
% perform differential analysis of metabolite data in drug pool vs all
% other pools at t=0 and t=12
fcMatrix_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
pMatrix_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
fcMatrix_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
pMatrix_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));

for j=1:length(ExperimentParameters.DrugNames)
    IdxPool = ExperimentParameters.PoolNumbers(ExperimentParameters.PoolingScheme(j,:)==1);
    IdxNOTPool = ExperimentParameters.PoolNumbers(ExperimentParameters.PoolingScheme(j,:)~=1);
    
    curData_t0 = totDataMat(totTime==0 & ismember(totDataPools, IdxPool),:);
    curData_t12 = totDataMat(totTime==12 & ismember(totDataPools, IdxPool),:);
    
    curCtrl_t0 = totDataMat(totTime==0 & ismember(totDataPools, IdxNOTPool),:);
    curCtrl_t12 = totDataMat(totTime==12 & ismember(totDataPools, IdxNOTPool),:);

    fcMatrix_t0(:,j) = (nanmean(curData_t0) ./...
                        nanmean(curCtrl_t0))';
    [~, pMatrix_t0(:,j)] = ttest2(curData_t0,...
                                  curCtrl_t0,...
                                  'VarType', 'equal');
                               
    fcMatrix_t12(:,j) = (nanmean(curData_t12) ./...
                            nanmean(curCtrl_t12))';
    [~, pMatrix_t12(:,j)] = ttest2(curData_t12,...
                                   curCtrl_t12,...
                                   'VarType', 'equal');
end
                                   
% adjust FDR
fdrMatrix_t0 = mafdr(pMatrix_t0(:), 'BHFDR', 1);
fdrMatrix_t0 = reshape(fdrMatrix_t0, size(fcMatrix_t0));
fdrMatrix_t12 = mafdr(pMatrix_t12(:), 'BHFDR', 1);
fdrMatrix_t12 = reshape(fdrMatrix_t12, size(fcMatrix_t12));
