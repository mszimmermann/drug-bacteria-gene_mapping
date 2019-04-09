function [fcMatrix_t0t12, fdrMatrix_t0t12, stdMatrix_t0t12,...
          rawIntensityMean_t0,rawIntensitySTD_t0,...
          rawIntensityMean_t12,rawIntensitySTD_t12] = workflow_DifferentialAnalysisT0T12(totDataMat,...
                                                                           ExperimentParameters,...
                                                                           totTime,...
                                                                           totDataPools)%,...
%                                                                            nobacteriaCTRLmat,...
%                                                                            nobacteriaCTRLtime,...
%                                                                            nobacteriaCTRLPools)
% perform differential analysis of metabolite data in drug pool vs all
% other pools at t=0 and t=12
pThreshold = 0.05;
fcThreshold = log2(2);

fcMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
pMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
stdMatrix_t0t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
fcMatrix_t0ctrl = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
pMatrix_t0ctrl = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
drugFC12to0_t0toCTRL = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));

rawIntensityMean_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensityMean_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensitySTD_t0 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));
rawIntensitySTD_t12 = zeros(size(totDataMat,2),length(ExperimentParameters.DrugNames));

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
    
    % controls are no-bacteria wells of the same pool      
%     curCtrl_t0 = nobacteriaCTRLmat(nobacteriaCTRLtime==0 & ismember(nobacteriaCTRLPools, IdxPool),:);
%   
%     fcMatrix_t0ctrl(:,j) = (nanmean(curData_t0) ./...
%                             nanmean(curCtrl_t0))';
%     [~, pMatrix_t0ctrl(:,j)] = ttest2(curData_t0,...
%                                    curCtrl_t0,...
%                                    'VarType', 'equal');
%     % calculate the difference between fold change t12 vs t0 
%     % and at t0 fold change between data vs ctrl
%     drugFC12to0_t0toCTRL(:,j) = (abs(log2(fcMatrix_t0t12(:,j))) - ...
%                                      abs(log2(fcMatrix_t0ctrl(:,j)))) .* ...
%                                      (log2(fcMatrix_t0t12(:,j))>0 & log2(fcMatrix_t0ctrl(:,j))>0);
end
% adjust FDR
fdrMatrix_t0t12 = mafdr(pMatrix_t0t12(:), 'BHFDR', 1);
fdrMatrix_t0t12 = reshape(fdrMatrix_t0t12, size(fcMatrix_t0t12));

%fdrMatrix_t0ctrl = mafdr(pMatrix_t0ctrl(:), 'BHFDR', 1);
%fdrMatrix_t0ctrl = reshape(fdrMatrix_t0ctrl, size(fcMatrix_t0ctrl));
 
% find drugs for which the difference between FC to ctrl at t0 and FC at
% t12 vs t0 is more than 2^5
% [~, fastDrugs_idx] = ind2sub(size(drugFC12to0_t0toCTRL), find(drugFC12to0_t0toCTRL<-5));
% fastDrugs_idx = unique(fastDrugs_idx);
% % for each of this drug, for each species choose max significant negative foldchange
% % between FC to ctrl at t0 and FC at t12 vs t0
% drugFC12to0_t0toCTRL_combined = log2(fcMatrix_t0t12);
% drugP12to0_t0toCTRL_combined = fdrMatrix_t0t12;
% 
% for idx=1:length(fastDrugs_idx)
%     j = fastDrugs_idx(idx);
%     change_fc = (log2(fcMatrix_t0ctrl(:,j))>=fcThreshold) &...
%                 (fdrMatrix_t0ctrl(:,j)<=pThreshold) &...
%                 (log2(fcMatrix_t0t12(:,j))<log2(fcMatrix_t0ctrl(:,j)));
%     drugFC12to0_t0toCTRL_combined(change_fc,j) = log2(fcMatrix_t0ctrl(change_fc,j));
%     drugP12to0_t0toCTRL_combined(change_fc,j) = fdrMatrix_t0ctrl(change_fc,j);
% end
% clear IdxPool curData_t0 curData_t12 curdrugidx                                  
% fcMatrix_t0t12 = drugFC12to0_t0toCTRL_combined;
% fdrMatrix_t0t12 = drugP12to0_t0toCTRL_combined;

