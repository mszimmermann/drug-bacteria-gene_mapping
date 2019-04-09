function [metFilter, metFragments] = filterDrugMetaboliteHits(metIntensity_t0,...
                                  metIntensity_t12,...
                                  metFoundFC_t0,...
                                  metFoundFC_t12,...
                                  drug_mets,...
                                  drugMass,...
                                  drugRT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                  
metFilter = zeros(size(drug_mets,1),5);
metFragments = zeros(size(drug_mets,1)); %indices of potential fragments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdBaseline = 5000;
thresholdFC = 1;
thresholdIntensity = 10000;
thresholdMZ = 0.002;
thresholdMassDefect = 0.2;
thresholdPPI = 0.02;
thresholdRT = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot mean intensity versus mean fold change over other plates
% #1 test intensity at t=0
metFound_Imean0 = metIntensity_t0;
metFound_Imean0(metFound_Imean0<=thresholdBaseline) = NaN;
metFound_FCmean0 = metFoundFC_t0;
metFound_FCmean0(metFound_FCmean0==0)=NaN;
metFound_FCmean0(isnan(metFound_Imean0))=NaN;

metFound_Imean0 = nanmean(metFound_Imean0,2);
metFound_FCmean0 = nanmean(log2(metFound_FCmean0),2);

metFilter(:,1) = metFound_Imean0>thresholdIntensity & metFound_FCmean0<thresholdFC;

% #2 test intensity at t=12
metFound_Imean12 = metIntensity_t12;
metFound_Imean12(metFound_Imean12<=thresholdBaseline) = NaN;
metFound_FCmean12 = metFoundFC_t12;
metFound_FCmean12(metFound_FCmean12==0)=NaN;
metFound_FCmean12(isnan(metFound_Imean12))=NaN;

metFound_Imean12 = nanmean(metFound_Imean12,2);
metFound_FCmean12 = nanmean(log2(metFound_FCmean12),2);

metFilter(:,2) = metFound_Imean12>thresholdIntensity & metFound_FCmean12<thresholdFC;

% #3 flag mass defects that are to far from the drug mass
drugidx = find((abs(drug_mets(:,2)-drugMass)<thresholdMZ) & ...
               (abs(drug_mets(:,3)-drugRT)<thresholdRT) );
if ~isempty(drugidx)
    if length(drugidx)>1
        [~, drugidx] = min(abs(drug_mets(drugidx,3)-drugRT));
    end
    massdefects = drug_mets(:,2) - floor(drug_mets(:,2));
    drug_massdefect = massdefects(drugidx);
    massdefects_diff = abs(massdefects-drug_massdefect);
    metFilter(:,3) = massdefects_diff>thresholdMassDefect;
else % drug is not detected untargetedly - use provided MZ and RT
    massdefects = drug_mets(:,2) - floor(drug_mets(:,2));
    drug_massdefect = drugMass - floor(drugMass);
    massdefects_diff = abs(massdefects-drug_massdefect);
    metFilter(:,3) = massdefects_diff>thresholdMassDefect;
end 
    
% #4 flag all metabolites that have different mass from drug,
% but the same RT
if ~isempty(drugidx)
    massdifference = abs(drug_mets(:,2) - drug_mets(drugidx,2));
    RTdifference = abs(drug_mets(:,3) - drug_mets(drugidx,3));
    metFilter(:,4) = massdifference > thresholdPPI*drug_mets(drugidx,2) &...
                     RTdifference < thresholdRT;
else % drug is not detected untargetedly - use provided MZ and RT
    massdifference = abs(drug_mets(:,2) - drugMass);
    RTdifference = abs(drug_mets(:,3) - drugRT);
    metFilter(:,4) = massdifference > thresholdPPI*drugMass &...
                     RTdifference < thresholdRT;
end
    
% #5 flag all masses that have been flagged at step #4 as masspec artifacts
% (if the drug has a real metabolite with a different RT, these metabolites 
% might be the massspec artifacts of drug metabolites
metFilter(:,5) = ismember(drug_mets(:,2), drug_mets(metFilter(:,4)==1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for each metabolite that has not been filtered, check whether there are
%potential fragments - metabolites with the same mass as drug fragments and
%the same RT as the metabolite
metidx = find(sum(metFilter,2)==0); %find metabolites that are not filtered out
for i=1:length(metidx)
    %check that the RT of the candidate is similar, 
    % and they occur in the sample sample at least once
    fragmentidx = metFilter(:,5) &...
                  abs(drug_mets(:,3)-drug_mets(metidx(i),3))<thresholdRT &...
                  (sum( (metIntensity_t0>thresholdIntensity) .*...
                    (repmat(metIntensity_t0(metidx(i),:),size(metIntensity_t0,1),1)>...
                     thresholdIntensity) ,2 ) |...
                   sum( (metIntensity_t12>thresholdIntensity) .*...
                    (repmat(metIntensity_t12(metidx(i),:),size(metIntensity_t12,1),1)>...
                     thresholdIntensity) ,2 ));
    metFragments(metidx(i),1:nnz(fragmentidx)) = find(fragmentidx);
end
metFragments(:, sum(metFragments,1)==0) = [];
end