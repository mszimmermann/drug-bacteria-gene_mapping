function reshapedMatrix = reshapePlateIntensities(linearData, sampleNames)
% reshape linear data as matrix of rows and columns
%extract rows and column numbers
rowNumbers = cellfun(@(x) str2double(x(strfind(x, '-R')+2:strfind(x, '-R')+3)), sampleNames);
columnNumbers = cellfun(@(x) str2double(x(strfind(x, '-C')+2:strfind(x, '-C')+3)), sampleNames);
rows_unique = unique(rowNumbers(~isnan(rowNumbers)));
cols_unique = unique(columnNumbers(~isnan(columnNumbers)));
reshapedMatrix = zeros(length(rows_unique), length(cols_unique));

[~, ~, cols_unique_indeces] = intersect(cols_unique, columnNumbers);
for i=1:length(rows_unique)
    % sum first row with all columns
    reshapedMatrix(i,:) = linearData(rowNumbers == rows_unique(i)) + ...
                            linearData(cols_unique_indeces)';
end