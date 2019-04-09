function [changingMets_merged,...
          changingMets_merged_idx,...
          changingMets_merged_idx_unique,...
          changingMets_merged_spectrum,...
          changingMets_merged_mass,...
          changingMets_merged_RT,...
          changingMets_merged_mass_delta,...
          changingMets_merged_RT_delta,...
          changingMets_merged_number] = merge_changing_metabolites(changingMets,...
                                                                   changingMetsSpectrum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge metabolites that are close in mass and RT
changingMets_merged = zeros(size(changingMets));
changingMets_mass = cellfun(@(x) str2double(x(1:strfind(x, '@')-1)), changingMets);
changingMets_RT = cellfun(@(x) str2double(x(strfind(x, '@')+1:end)), changingMets);
massThreshold = 0.002;
massThresholdSpectrumCoeff=5; %for composite spectra, use relaxed massThresholdSpectrumCoeff*massThreshold
RTthreshold = 0.2;
for i=1:length(changingMets_merged)
    if changingMets_merged(i)==0
        % find candidates with similar mass and RT
        sameidx = find( (abs(changingMets_mass-changingMets_mass(i))<=massThreshold) &...
                    (abs(changingMets_RT-changingMets_RT(i))<=RTthreshold) );
        % expand the candidates by adding candidates of smalles and larges mz from the list
        [max_mass, idx] = max(changingMets_mass(sameidx));
        max_rt = changingMets_RT(sameidx(idx));
        [min_mass, idx] = min(changingMets_mass(sameidx));
        min_rt = changingMets_RT(sameidx(idx));
        sameidx_max = find( (abs(changingMets_mass-max_mass)<=massThreshold) &...
                    (abs(changingMets_RT-max_rt)<=RTthreshold) );
        sameidx_min = find( (abs(changingMets_mass-min_mass)<=massThreshold) &...
                    (abs(changingMets_RT-min_rt)<=RTthreshold) );        
        sameidx = unique([sameidx; sameidx_min; sameidx_max]);
        if length(sameidx)>1
            curSpectra = changingMetsSpectrum(sameidx);
            % find ions that have more than one overlap with the i-th
            baseSpectrum = strsplit(curSpectra{sameidx==i}, {'(',',',')'});
            baseSpectrum = cellfun(@(x) str2double(x), baseSpectrum(2:2:end));
            baseSpectrum(isnan(baseSpectrum)) = [];
            % for each same ion suspect, calculate how many ions from spectrum
            % overlap
            spectrumOverlap = zeros(length(curSpectra),1);
            for j=1:length(curSpectra)
                curSpectraMZ = strsplit(curSpectra{j}, {'(',',',')'});
                curSpectraMZ = cellfun(@(x) str2double(x), curSpectraMZ(2:2:end));
                curSpectraMZ(isnan(curSpectraMZ)) = [];
                curSpectra{j} = curSpectraMZ;
                mzdiff = bsxfun(@(x,y) x-y, baseSpectrum', curSpectra{j});
                spectrumOverlap(j) = nnz(sum(abs(mzdiff)<=massThresholdSpectrumCoeff*massThreshold,2));
            end
            sameidx = sameidx(spectrumOverlap>1);
        end
        % add the ion to existing group if it already exists
        if nnz(changingMets_merged(sameidx))
            existing_idx = changingMets_merged(sameidx);
            existing_idx = unique(existing_idx(existing_idx~=0));
            if length(existing_idx)==1
                changingMets_merged(sameidx) = unique(existing_idx);
            else
                % lump all ions that match the current one together
                for j=1:length(existing_idx)
                    changingMets_merged(changingMets_merged==existing_idx(j)) = i;
                end
                changingMets_merged(sameidx) = i;
            end
        else
            changingMets_merged(sameidx) = i;
        end
    end
end
changingMets_merged_idx = changingMets_merged;
changingMets_merged_idx_unique = unique(changingMets_merged_idx);
% calculate average mass and RT of all merged metabolites
changingMets_merged = cell(size(changingMets_merged_idx_unique));
changingMets_merged_spectrum = cell(size(changingMets_merged_idx_unique));
changingMets_merged_mass = zeros(size(changingMets_merged_idx_unique));
changingMets_merged_RT = zeros(size(changingMets_merged_idx_unique));
changingMets_merged_mass_delta = zeros(size(changingMets_merged_idx_unique));
changingMets_merged_RT_delta = zeros(size(changingMets_merged_idx_unique));
changingMets_merged_number = zeros(size(changingMets_merged_idx_unique));
for i=1:length(changingMets_merged_idx_unique)
    % merge the ions and take the median mass as the reference
    merged_idx = find(changingMets_merged_idx==changingMets_merged_idx_unique(i));
    changingMets_merged_mass(i) = median(changingMets_mass(merged_idx));
    [~, med_idx] = min(abs(changingMets_mass(merged_idx)-changingMets_merged_mass(i)));
    % get the median ion's RT, calculate mass deltas and RT deltas
    changingMets_merged_RT(i) = changingMets_RT(merged_idx(med_idx));
    changingMets_merged_mass_delta(i) = max(abs(changingMets_mass(merged_idx) -...
                                                changingMets_merged_mass(i)));
    changingMets_merged_RT_delta(i) = max(abs(changingMets_RT(merged_idx)-...
                                                changingMets_merged_RT(i)));                                            
    changingMets_merged{i} = [num2str(changingMets_merged_mass(i)), '@',...
                              num2str(changingMets_merged_RT(i))];
                          
    changingMets_merged_spectrum{i} = changingMetsSpectrum{merged_idx(med_idx)};
    changingMets_merged_number(i) = length(merged_idx);
end

clear sameidx spectrumOverlap mzdiff curSpectra curSpectraMZ
clear baseSpectrum max_mass idx max_rt min_mass min_rt sameidx_max sameidx_min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
