function [maxSlope, maxSlopeB] = calculateConversionSlope(curTime_unique, meanPerTime, timeWindow, SlopeDir)

    curCoefs = zeros(length(curTime_unique)-timeWindow,2);
%     figure;
%     plot(curTime_unique, log2(meanPerTime));
%     hold on
   
    for i=1:length(curTime_unique)-timeWindow
        timeFitThreshold = curTime_unique(i):curTime_unique(i+timeWindow);
        curfit = polyfit(curTime_unique(ismember(curTime_unique,timeFitThreshold)),...
                         meanPerTime(ismember(curTime_unique,timeFitThreshold)), 1);
        curCoefs(i,1) = curfit(1);
        curCoefs(i,2) = curfit(2);

%         plot(curTime_unique(ismember(curTime_unique,timeFitThreshold)),...
%                  curfit(1)*curTime_unique(ismember(curTime_unique,timeFitThreshold))+curfit(2), '--')

    end
    if SlopeDir>0
        [~, idx] = max(curCoefs(:,1));
    else
        [~, idx] = min(curCoefs(:,1));
    end
    maxSlope = curCoefs(idx,1);
    maxSlopeB = curCoefs(idx,2);
%    timeFitThreshold = curTime_unique(idx):curTime_unique(idx+timeWindow);
%     plot(curTime_unique(ismember(curTime_unique,timeFitThreshold)),...
%                  maxSlope*curTime_unique(ismember(curTime_unique,timeFitThreshold))+maxSlopeB, '--', 'LineWidth',2)

       