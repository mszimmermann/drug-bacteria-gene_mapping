function [bestP,bestStat,bestFC,bestPident] = metabolizer_gene_enrichment(...
                                                        curAdjustedFC_pident,...
                                                        curAdjustedFC,...
                                                        sortdir)
if size(curAdjustedFC,1) == 1
    curAdjustedFC = curAdjustedFC';
end
if size(curAdjustedFC_pident,1) == 1
    curAdjustedFC_pident = curAdjustedFC_pident';
end

[~, sort_pident] = sort(curAdjustedFC_pident, 'descend');

sort_pident = sort_pident(curAdjustedFC_pident(sort_pident)>0);

% leave only unique values for FC and pident
sort_pident_val = curAdjustedFC_pident(sort_pident);
[sort_pident_val, idx] = unique(sort_pident_val, 'stable');
sort_pident = sort_pident(idx);
% include pident thresholds down to 50%
idx = (sort_pident_val>=50);
sort_pident_val = sort_pident_val(idx);
sort_pident = sort_pident(idx);
%%%%%%%%%%%%% FC threshold are not defined by the data but by the user
sort_fc_val = fliplr(0.2:0.2:0.8); % fold change values 20% 40% 60% 80%
sort_fc = 1:length(sort_fc_val);

nSpecies = nnz(~isnan(curAdjustedFC_pident));
enrichmentPmatrix = ones(length(sort_fc), length(sort_pident));
enrichmentStatmatrix = cell(length(sort_fc), length(sort_pident));
for pi = 1:length(sort_pident)
    for fci = 1:length(sort_fc)
        nGroupCountainingGene = nnz(curAdjustedFC_pident>=sort_pident_val(pi));% this is K - number of species containing gene
        
        if sortdir==-1
            % this is k - number of species containing the gene which are also metabolizers
            nChangingGroup = nnz(curAdjustedFC_pident>=sort_pident_val(pi) & ...
                                 curAdjustedFC<=sort_fc_val(fci));
            % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
            % where N - total number of genes, n - size of the group
            nChanging = nnz(curAdjustedFC<=sort_fc_val(fci));
        else
            % this is k - number of species containing the gene which are also metabolizers
            nChangingGroup = nnz(curAdjustedFC_pident>=sort_pident_val(pi) & ...
                                 curAdjustedFC>=sort_fc_val(fci));
            % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
            % where N - total number of genes, n - size of the group
            nChanging = nnz(curAdjustedFC>=sort_fc_val(fci));
        end

        score = sum( hygepdf(nChangingGroup:nGroupCountainingGene,...
                             nSpecies,...
                             nGroupCountainingGene,...
                             nChanging));
        enrichmentPmatrix(fci,pi) = score;
        enrichmentStatmatrix{fci,pi} = [nChangingGroup,nSpecies, nGroupCountainingGene,nChanging];
    end
end
[best_fc,best_pid] = find(enrichmentPmatrix==min(min(enrichmentPmatrix)));
if length(best_fc) > 1
    best_fc = best_fc(1);
    best_pid = best_pid(1);
end
bestP = enrichmentPmatrix(best_fc, best_pid);
bestStat = enrichmentStatmatrix{best_fc, best_pid};
bestFC = sort_fc_val(best_fc);
bestPident = sort_pident_val(best_pid);


