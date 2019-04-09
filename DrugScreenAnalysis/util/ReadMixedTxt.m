function TempStruct = ReadMixedTxt(datafile, delimiter)
  TempStruct = struct;
  fid = fopen(datafile,'r');   %# Open the file
  lineArray = cell(50000,1);     %# Preallocate a cell array (ideally slightly   %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
      lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
  end
  %replace empty cells by 1s
  lineArray(cellfun(@isempty,lineArray))={'1'};
  % parse data into structure
  % import raw data
  Idx = cellfun(@(x) ~isempty(strfind(x, 'Area')), lineArray(2,:));
  TempStruct.IntensitiesRaw =  cellfun(@(x) str2num(x), lineArray(3:end, Idx));
  TempStruct.IntensitiesRaw(TempStruct.IntensitiesRaw==0) = 1;
  % import RT data
  Idx = cellfun(@(x) ~isempty(strfind(x, 'RT')), lineArray(2,:));
  TempStruct.RT =  cellfun(@(x) str2num(x), lineArray(3:end, Idx));
  % import compounds
  TempStruct.Compounds =  lineArray(1, Idx);
  TempStruct.SampleNames =  lineArray(3:end, 4);
  
  % import sample names
  
  

%    % Import composite spectrum
%    CompositeSpectrum = cellfun(@(x) ~isempty(strfind(x, 'CompositeSpectrum')), lineArray(5, :));
%    TempStruct.CompositeSpectrum =  lineArray(5, CompositeSpectrum);
%    % Import RT
%    RT = cellfun(@(x) ~isempty(strfind(x, 'Retention')), lineArray(5, :));
%    TempStruct.RT =  cellfun(@(x) str2num(x),lineArray(6:end, RT));
%    % Import Frequency
%    Frequency = cellfun(@(x) ~isempty(strfind(x, 'Frequency')), lineArray(5, :));
%    TempStruct.Frequency =  cellfun(@(x) str2num(x),lineArray(6:end, Frequency));
%    % Import Masses
%    Masses = cellfun(@(x) ~isempty(strfind(x, 'Mass')), lineArray(5, :));
%    TempStruct.Masses =  cellfun(@(x) str2num(x),lineArray(6:end, Masses));
  
    
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    % find outlayer data due to bad injections, if sum of all compounds is
%    % lower than the mean of all 48 data sets minus three times their std (=
%    % threshold value => search for injection failures
%    Threshold = mean(sum(TempStruct.IntensitiesRaw)-(3*(std(sum(TempStruct.IntensitiesRaw)))));
%    TempStruct.ExcludedDatasets = find(sum(TempStruct.IntensitiesRaw) < Threshold);
%    TempStruct.IntensitiesNorm = TempStruct.IntensitiesRaw;
%    TempStruct.IntensitiesNorm(:,TempStruct.ExcludedDatasets) = [];
%    TempStruct.Injections(TempStruct.ExcludedDatasets) = [];
%    % TempStruct.IntensitiesNorm = quantilenorm(TempStruct.IntensitiesNorm);%
%    % quantile normalization over all the time points
%    
%     
%    % exclude ions that were unique to the excuded dataset in all
%    % parameters
%    TempStruct.IdxExcluded = find(sum(TempStruct.IntensitiesNorm,2) ==0);
%    TempStruct.Compounds(TempStruct.IdxExcluded,:) =  [];
%    TempStruct.CompositeSpectrum(TempStruct.IdxExcluded,:) = [];
%    TempStruct.RT(TempStruct.IdxExcluded,:) = [];
%    TempStruct.Frequency(TempStruct.IdxExcluded,:) =  [];
%    TempStruct.Masses(TempStruct.IdxExcluded,:) = [];
%    
%    % find missing and excluded pools per time point
%    TempStruct.TimeIdx = zeros(size(ExperimentParameters.TimePoints,2), size(TempStruct.Injections,2)); % Idx of injections of same time point
%    TempStruct.PoolIncluded = zeros(size(ExperimentParameters.TimePoints,2), size(ExperimentParameters.PoolNo,2)); % pools per time point
%    
%    for i = 1: size(ExperimentParameters.TimePoints,2)
%      TempStruct.TimeIdx(i,:) = cell2mat(cellfun(@(x) ~isempty(strfind(x, ExperimentParameters.TimePoints{i})), TempStruct.Injections, 'UniformOutput', false));
%      TempInjections = TempStruct.Injections(TempStruct.TimeIdx(i,:)==1);
%      TempPools = cellfun(@(x) x(12:13), TempInjections, 'UniformOutput', false);
%      TempPools = cellfun(@(x) str2num(x), TempPools);
%      TempStruct.PoolIncluded(i,:) = ismember(ExperimentParameters.PoolNo, TempPools);
%    end  
%    
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  % Normalize intensities by ions that are measured in all samples on one community  
%   IdxAll = find(TempStruct.Frequency == size(TempStruct.IntensitiesNorm,2));
%   DataTemp = TempStruct.IntensitiesNorm(IdxAll, :);
%   
%   FC_IS = zeros(length(IdxAll), size(TempStruct.IntensitiesNorm, 2));
%   MeanAll = mean(DataTemp,2);
%   
%   for i = 1:length(MeanAll)
%      FC_IS(i,:) = TempStruct.IntensitiesNorm(IdxAll(i),:) ./ MeanAll(i);
%   end
%    
%   MeanFC_IS = median(FC_IS(1:2,:),1);
%    
%    for i = 1: length(MeanFC_IS)
%        Datatemp = TempStruct.IntensitiesNorm(:,i);
%        for p = 1: length(Datatemp)
%            if Datatemp(p) > 1
%            Datatemp(p) =  Datatemp(p)/MeanFC_IS(i);
%            else
%            end
%        end
%    end
%   
%   
%    
%    
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  %   %perform quantile normalization per time point
%  %   for t = 1:length(ExperimentParameters.TimePoints)
%  %       Datatemp = TempStruct.IntensitiesNorm(:,(TempStruct.TimeIdx(t,:)==1));
%  %       Datatemp = quantilenorm(Datatemp);
%  %       TempStruct.IntensitiesNorm(:,(TempStruct.TimeIdx(t,:)==1)) = Datatemp;
%  %   end
%  
%  
%   
%    StdIdx = zeros(size(InternalStd.StdRT,1), 1);
%    % find internal STD for normalization 
%    for i = 1:size(InternalStd.StdRT,1)
%        RTWindow = find(TempStruct.RT > (InternalStd.StdRT(i)-0.2) & TempStruct.RT < (InternalStd.StdRT(i)+0.2));
%        MinMass = find(abs(TempStruct.Masses(RTWindow)- InternalStd.StdMasses(i)) == min(abs((TempStruct.Masses(RTWindow) - InternalStd.StdMasses(i)))));
%        StdIdx(i) = RTWindow(MinMass);  
%    end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % normalize by Internal STD
%   FC_IS = zeros(length(StdIdx), size(TempStruct.IntensitiesNorm, 2));
%   MeanIS = mean(TempStruct.IntensitiesNorm(StdIdx,:),2);
%   for i = 1:length(MeanIS)
%   FC_IS(i,:) = TempStruct.IntensitiesNorm(StdIdx(i),:) ./ MeanIS(i);
%   end
%   
%   MeanFC_IS = mean(FC_IS(1:2,:),1);
%   
%   for i = 1: length(MeanFC_IS)
%       Datatemp = TempStruct.IntensitiesNorm(:,i);
%       for p = 1: length(Datatemp)
%           if Datatemp(p) > 1
%           Datatemp(p) =  Datatemp(p)/MeanFC_IS(i);
%           else
%           end
%       end
%   end
      

% figure      
% for p = 1:5
% subplot (2,3,p)
% plot(1:216, TempStruct.IntensitiesNorm(StdIdx(p),:))
% end
%   
  
  
end