function [maxMZ] = parseMaxCompositeSpectra(curspectra)
% return mass and RT of the largest ion intensity
curIons = strsplit(curspectra{1},{',' ')'});
curIons = cellfun(@(x) strrep(x, '(', ''), curIons, 'unif', 0);
curIons = cellfun(@(x) strtrim(x), curIons, 'unif', 0);
curIons(cellfun(@(x) isempty(x), curIons)) = [];
curIons = cellfun(@(x) str2double(x), curIons);
[~, maxidx] = max(curIons(2:2:end));
maxMZ = curIons(maxidx*2-1);
