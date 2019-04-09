function outliers = detect_one_outlier(curData, outlierSTDcutoff, outlierMEANcutoff)
% function that detects one outlier out of the vector of values 
% based on the difference of the mean when the outlier is removed
% and the diffeence between uotlier and the mean of other replicates 
% in terms of STD
curData_std_diff = zeros(size(curData));
curData_mean_diff = zeros(size(curData));
for k=1:length(curData)
    curData_std_diff(k) = abs(mean(curData(setdiff(1:length(curData),k),:))-curData(k,:))/std(curData(setdiff(1:length(curData),k),:));
    curData_mean_diff(k) = abs(mean(curData(setdiff(1:length(curData),k),:))-mean(curData))/mean(curData);
end
outliers = curData_std_diff>outlierSTDcutoff & curData_mean_diff>outlierMEANcutoff;
%if there are more than 1 outlier, do not call the replicates outliers
if nnz(outliers)>1
    outliers = zeros(size(curData));
end
       