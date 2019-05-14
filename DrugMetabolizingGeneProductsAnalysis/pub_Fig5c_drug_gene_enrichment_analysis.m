%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig5c_drugmet_gene_enrichment_analysis
% read drug metabolism data and gene-genomes identity BLAST file
% perform gene enrichment analysis among drug-metabolizing species with
% pident to the identified drug.metabolizing gene products (with sliding
% thresholds both for drug consumption and gene pident)
% plot enrichment heatmaps (for each gene and combination of genes)
% and save them to pdf files, 
% save enrichment results to table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% will use the following variables:
% infileTable14BLAST % BLAST table of identified genes and genomes of single species
% infile_gene_drug_matrix % Table of curated gene-drug interactions
% outfile_table_FCclustergram % metabolized drug fol changes in clustergram order
% infile_single_species % single species information
% outfile_fig5c_drug_gene_enrichment
% outfile_fig5c_combined_enrcihment
% outfile_table15_gene_enrcihment_table

% read blast information of genes and genomes
genome_blast = readtable(infileTable14BLAST);           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read gene_drug_metabolite table
gene_drug_table = readtable(infile_gene_drug_matrix);
% extract information on gene_drugs (find parent columns in the matrix)
column_types = varfun(@class,gene_drug_table,'OutputFormat','cell');
column_parent = zeros(size(column_types));
for i=1:length(column_types)
    if isequal(column_types{i}, 'cell')
        column_parent(i) = contains(lower(gene_drug_table{1,i}),'parent');
    end
end
column_metabolite = column_parent==0;
column_metabolite(ismember(gene_drug_table.Properties.VariableNames, 'Gene') |...
                  cellfun(@(x) contains(x, 'Var'), gene_drug_table.Properties.VariableNames)) = 0;
col_met = find(column_metabolite);
for i=1:length(col_met)
    if ~isnumeric(gene_drug_table.(gene_drug_table.Properties.VariableNames{col_met(i)}))
        gene_drug_table.(gene_drug_table.Properties.VariableNames{col_met(i)}) = ...
            cellfun(@(x) str2double(x), gene_drug_table.(gene_drug_table.Properties.VariableNames{col_met(i)}));
    end
end
        
% create list of gene-drug pairs
gene_drug_list = cell(size(gene_drug_table,1)*size(gene_drug_table,2),2);
idx = 1;
for i=2:size(gene_drug_table,1)
    curdrugmat = cellfun(@(x) str2double(x), gene_drug_table{i,column_parent==1});
    curdrugs = gene_drug_table.Properties.VariableNames(column_parent==1);
    curdrugs = curdrugs(curdrugmat==1);
    curmetsmat = gene_drug_table{i,column_metabolite==1};
    curmets = cellfun(@(x) x(1:strfind(x, '_')-1), gene_drug_table.Properties.VariableNames(column_metabolite==1), 'unif',0);
    curmets = unique(curmets(curmetsmat==1));
    curdrugs = union(curdrugs, curmets);
    gene_drug_list(idx:idx+length(curdrugs)-1,:) = [repmat(gene_drug_table.Gene(i),size(curdrugs)); curdrugs]';
    idx = idx+length(curdrugs);
end
gene_drug_list(idx:end,:) = [];
% add space before acetate that got lost during variable import
gene_drug_list(:,2) = cellfun(@(x) strrep(x, 'Acetate', ' acetate'), gene_drug_list(:,2), 'unif', 0);
gene_drug_list(:,2) = cellfun(@(x) strrep(x, 'Diacetate', ' diacetate'), gene_drug_list(:,2), 'unif', 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load clustergram of changing drugs to get drugs and species
drugClustergram = readtable(outfile_table_FCclustergram,...
                            'HeaderLines', 1);
drugAnalysisColumnNames = drugClustergram.Properties.VariableNames;
FCcolumns = cellfun(@(x) contains(x, 'FC') &...
                         ~contains(x, 'FCSTD') &...
                         ~contains(x, 'DrugAdaptive'),...
                         drugAnalysisColumnNames);
FCSTDcolumns = cellfun(@(x) contains(x, 'FCSTD'),...
                         drugAnalysisColumnNames);
PFDRcolumns = cellfun(@(x) contains(x, 'p_FDR'),...
                         drugAnalysisColumnNames);                     
drugFC12to0_t0toCTRL_combined = drugClustergram{:,FCcolumns}';
drugFC12to0_t0toCTRL_STDcombined = drugClustergram{:,FCSTDcolumns}';
drugP12to0_t0toCTRL_combined = drugClustergram{:,PFDRcolumns}';
drugFC12to0_Species = drugAnalysisColumnNames(FCcolumns);
% replace FC with empty
drugFC12to0_Species = cellfun(@(x) strrep(x, 'FC', ''), drugFC12to0_Species, 'unif',0);
% raformat strain names
for si = 1:length(drugFC12to0_Species)
    curspace = find(isstrprop(drugFC12to0_Species{si},'upper') | isstrprop(drugFC12to0_Species{si},'digit'));
    if length(curspace)>2
        drugFC12to0_Species{si} = [drugFC12to0_Species{si}(1:curspace(2)-1)...
                          ' '...
                          lower(drugFC12to0_Species{si}(curspace(2):curspace(3)-1))...
                          ' '...
                          drugFC12to0_Species{si}(curspace(3):end)];
    elseif length(curspace)>1
        drugFC12to0_Species{si} = [drugFC12to0_Species{si}(1:curspace(2)-1)...
                          ' '...
                          lower(drugFC12to0_Species{si}(curspace(2):end))];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load single spesies information 
AssaySpecies = table2cell(readtable(infile_single_species));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qseqidx = cellfun(@(x) contains(x, 'qseqid'), genome_blast.Properties.VariableNames);
sseqidx = cellfun(@(x) contains(x, 'sseqid'), genome_blast.Properties.VariableNames);
pidentidx = cellfun(@(x) contains(x, 'pident'), genome_blast.Properties.VariableNames);

genome_name_strain = genome_blast.(genome_blast.Properties.VariableNames{sseqidx});
genome_name_strain = cellfun(@(x) strrep(x,'ATCC ', 'ATCC'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'DSM ', 'DSM'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'CCUG ', 'CCUG'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'nexilis', 'nexile'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'reuteri CF48-3A', 'reuteri CF48-3A BEI HM-102'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'xylanisolvens CL03T12C04', 'xylanisolvens CL03T12C04 DSM18836'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'indistinctus YIT 12060', 'indistinctus YIT 12060 DSM 22520'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'biformis DSM3989', 'biforme DSM3989'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'indistinctus YIT 12060', 'indistinctus YIT 12060 DSM22520'),genome_name_strain,'unif',0);
genome_name_strain = cellfun(@(x) strrep(x,'vadensis strain DSM14823', 'vadensis strain DSM 14823 ATCC BAA-548'),genome_name_strain,'unif',0);

% match assay species with NP genome IDs
AssaySpecies_matched = zeros(size(AssaySpecies,1),10); 
for i=1:size(AssaySpecies,1)
    if ~isempty(AssaySpecies{i,5})
        if isnumeric(AssaySpecies{i,5})
            curidx = find(cellfun(@(x) contains(lower(x), lower(num2str(AssaySpecies{i,5}))) &...
                                       contains(lower(x), lower(AssaySpecies{i,4})),...
                                       genome_name_strain));
        else
            curidx = find(cellfun(@(x) contains(lower(x), lower(AssaySpecies{i,5})) &...
                                       contains(lower(x), lower(AssaySpecies{i,4})),...
                                       genome_name_strain));
        end
        AssaySpecies_matched(i, 1:length(curidx)) = curidx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% match drug clustergram species with S numbers
AssaySpecies(:,5) = cellfun(@(x) strrep(x, ' ', ''), AssaySpecies(:,5), 'unif', 0);
AssaySpecies(:,5) = cellfun(@(x) strrep(x, '-', ''), AssaySpecies(:,5), 'unif', 0);
AssaySpecies(:,5) = cellfun(@(x) strrep(x, '_', ''), AssaySpecies(:,5), 'unif', 0);
AssaySpecies(:,5) = cellfun(@(x) strrep(x, '(', ''), AssaySpecies(:,5), 'unif', 0);
AssaySpecies(:,5) = cellfun(@(x) strrep(x, ')', ''), AssaySpecies(:,5), 'unif', 0);
AssaySpecies(:,4) = cellfun(@(x) lower(x), AssaySpecies(:,4), 'unif', 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drugFC12to0_Species = cellfun(@(x) strrep(x, '_', ''), drugFC12to0_Species, 'unif', 0);
drugFC12to0_SpeciesS = cell(size(drugFC12to0_Species));
for i=1:length(drugFC12to0_Species)
    curS = AssaySpecies( cellfun(@(x) contains(drugFC12to0_Species{i}, x), AssaySpecies(:,4)) &...
                         cellfun(@(x) contains(drugFC12to0_Species{i}, x), AssaySpecies(:,5)), 1);
    if ~isempty(curS)
        drugFC12to0_SpeciesS{i} = sprintf('S%03d', curS{1});
    else
        drugFC12to0_SpeciesS{i} = 'S000';
    end
end
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'longum'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'longum'),AssaySpecies(:,4)), 1});
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'Anaerostipes sp'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'sp.'), AssaySpecies(:,4)) &...
                                  cellfun(@(x) contains(x, 'Anaerostipes'), AssaySpecies(:,3)), 1});
drugFC12to0_SpeciesS{cellfun(@(x) isequal(x, 'Clostridium sp'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'sp.'), AssaySpecies(:,4)) &...
                                  cellfun(@(x) contains(x, 'Clostridium'), AssaySpecies(:,3)), 1});
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'reuteri'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'reuteri'),AssaySpecies(:,4)), 1});
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'WH2'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'wh2'),AssaySpecies(:,4)), 1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each genome get the corresponding species
genome_id_species_number = zeros(size(genome_name_strain,1),1);
for i=1:size(genome_name_strain,1)
    [rowi, coli] = find(AssaySpecies_matched==i);
    if rowi
        curS = cell2mat(AssaySpecies(rowi,1));
        genome_id_species_number(i) = curS;
    end
end
genome_id_species_S = arrayfun(@(x) sprintf('S%03d',x), genome_id_species_number, 'unif', 0);

species_no_genomes = cellfun(@(x) ~ismember(x, genome_id_species_S), drugFC12to0_SpeciesS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each drug and gene pair perform enrichment analysis
% with sliding threshold for drug fold change
% AND sliding threshold for pident

% get names of drugs for which gene is identified
drug_names_unique = gene_drug_list(:,2);
% remove underscores if there are any in drug names
for i=1:length(drug_names_unique)
    if contains(drug_names_unique{i}, '_')
        drug_names_unique{i} = drug_names_unique{i}(1:strfind(drug_names_unique{i},'_')-1);
    end
end   
drug_names_unique = unique(drug_names_unique);

metGenes = unique(gene_drug_list(:,1));%unique(genome_blast.qseqid);
metGenes(cellfun(@(x) ~contains(x, '_'), metGenes)) = [];
        
metGenesDrugs_enrP = ones(length(metGenes), length(drug_names_unique));
metGenesDrugs_enrStat = cell(length(metGenes), length(drug_names_unique));
metGenesDrugs_enrFC = zeros(length(metGenes), length(drug_names_unique));
metGenesDrugs_enrPident = zeros(length(metGenes), length(drug_names_unique));
% combined statistics for max identity to either of the metabolizing genes
metGenesDrugs_combined_enrP = ones(length(drug_names_unique),1);
metGenesDrugs_combined_enrStat = cell(length(drug_names_unique),1);
metGenesDrugs_combined_enrFC = zeros(length(drug_names_unique),1);
metGenesDrugs_combined_enrPident = zeros(length(drug_names_unique),1);
for drug_i =1:length(drug_names_unique)
    %find genes that metabolize the drug
    curgeneidx = ismember(gene_drug_list(:,2),drug_names_unique(drug_i));
    
    % combined pident for all drug genes
    drug_clustidx = ismember(lower(drugClustergram.DrugName), lower(drug_names_unique(drug_i)));
    curAdjustedFC = drugFC12to0_t0toCTRL_combined(:, drug_clustidx);
    curAdjustedFC = 1-2.^curAdjustedFC ;
    
    if nnz(curgeneidx)
        curgenes = unique(gene_drug_list(curgeneidx,1));
        combinedGenePident = zeros(length(curAdjustedFC),length(curgenes));  

        for gene_i = 1:length(curgenes)
            curgenomesidx = cellfun(@(x) contains(x, curgenes(gene_i)),...
                                          genome_blast.(genome_blast.Properties.VariableNames{qseqidx}));
            if nnz(curgenomesidx)
                curPident = genome_blast.(genome_blast.Properties.VariableNames{pidentidx});
                curPident = curPident(curgenomesidx);

                curgenomes_AssayS = genome_id_species_S(curgenomesidx);
                % leave only one pident per strain (max value)
                curgenomes_AssayS_unique = unique(curgenomes_AssayS);
                curgenomes_AssayS_unique(ismember(curgenomes_AssayS_unique, 'S000'))=[];
                for assay_i = 1:length(curgenomes_AssayS_unique)
                    curidx = find(ismember(curgenomes_AssayS,curgenomes_AssayS_unique{assay_i}));
                    curpident = curPident(curidx);
                    [~, maxidx] = max(curpident);
                    curidx(maxidx) = [];
                    curgenomes_AssayS(curidx) = {'S000'};
                end
                % match species with drug FC
                curAdjustedFC_pident = cellfun(@(x) curPident(ismember(curgenomes_AssayS, x)),...
                                                    drugFC12to0_SpeciesS, 'unif', 0);
                curAdjustedFC_pident(cellfun(@(x) isempty(x), curAdjustedFC_pident)) = {0};
                curAdjustedFC_pident = cell2mat(curAdjustedFC_pident);
                curAdjustedFC_pident(species_no_genomes) = nan;

                combinedGenePident(:, gene_i) = curAdjustedFC_pident;
                % perform gene nerichment with sliding thresholds
                [bestP,bestStat,bestFC,bestPident] = metabolizer_gene_enrichment(...
                                                            curAdjustedFC_pident(~isnan(curAdjustedFC_pident)),...
                                                            curAdjustedFC(~isnan(curAdjustedFC_pident)), 1);

                metGene_idx = ismember(metGenes, curgenes{gene_i});
                metGenesDrugs_enrP(metGene_idx,drug_i) = bestP;
                metGenesDrugs_enrStat{metGene_idx,drug_i} = bestStat;
                metGenesDrugs_enrFC(metGene_idx,drug_i) = bestFC;
                metGenesDrugs_enrPident(metGene_idx,drug_i) = bestPident;
            end
        end
        combinedGenePident = max(combinedGenePident,[],2);
        % perform gene nerichment with sliding thresholds
        % for combined pident (max of all genes metabolizing the drug)
        [bestP,bestStat,bestFC,bestPident] = metabolizer_gene_enrichment(...
                                                    combinedGenePident(~isnan(curAdjustedFC_pident)),...
                                                    curAdjustedFC(~isnan(curAdjustedFC_pident)), 1);

        metGenesDrugs_combined_enrP(drug_i) = bestP;
        metGenesDrugs_combined_enrStat{drug_i} = bestStat;
        metGenesDrugs_combined_enrFC(drug_i) = bestFC;
        metGenesDrugs_combined_enrPident(drug_i) = bestPident;
    end
    fprintf('Done with drug %s\n', drug_names_unique{drug_i})
end

% perfrom FDR correction for multiple hypotheses testing across all p-values
metGenesDrugs_enrFDR = mafdr(metGenesDrugs_enrP(:), 'BHFDR', 1);
metGenesDrugs_enrFDR = reshape(metGenesDrugs_enrFDR, size(metGenesDrugs_enrP));

imageFDR = metGenesDrugs_enrFDR;
selectCols = 1:size(imageFDR,2);
imageFDR(imageFDR>0.05) = 1;
selectRows = 1:size(imageFDR,1);
imageFDR = imageFDR(selectRows, selectCols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display clustergram of the enrichment fdr values
cmaplength = 6;
clustdist = 'euclidean';
mycmap = repmat(linspace(0, 1,cmaplength)',1,3);

set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
cluster_matrix = -log10(imageFDR);
cmaprange = 6;

cgo=clustergram(cluster_matrix,...
                 'DisplayRange', cmaprange,...
                 'RowLabels', metGenes(selectRows),...
                 'ColumnLabels', drug_names_unique(selectCols),...
                 'Symmetric', 0,...
                 'ColumnPdist', clustdist, 'RowPdist', clustdist,...
                 'Colormap', mycmap);

set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 6)
 
% insert colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
cb  = findobj(gcf,'Tag','HeatMapColorbar');
cb.Label.String = 'Enrichment p(FDR)';

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    outfile_fig5c_drug_gene_enrichment)

                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metGenesDrugs_combined_enrFDR = mafdr(metGenesDrugs_combined_enrP, 'BHFDR', 1);
% plot enrichment for all genes per drug in the same order as clustergram
imageFDR = metGenesDrugs_combined_enrFDR(selectCols);
imageFDR(imageFDR>0.05)=1;
colLabels = drug_names_unique(selectCols);
[~, ~, idx] = intersect(cgo.ColumnLabels, colLabels, 'stable');
imageFDR = imageFDR(idx)';
colLabels = colLabels(idx);

fig = figure;
imagesc(-log10(imageFDR))
set(gca, 'XTick', 1:nnz(selectCols))
set(gca, 'XTickLabel', colLabels)
xticklabel_rotate([], 90)
colormap(mycmap)
caxis([0 cmaprange]);
colorbar;
cb  = findobj(gcf,'Tag','Colorbar');
cb.Label.String = 'Enrichment p(FDR)';

title('Combination of genes enrichment among drug metabolizers')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
        outfile_fig5c_combined_enrcihment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(outfile_table15_gene_enrcihment_table,'w');
printDrugs = 1:length(drug_names_unique);
fprintf(fid,'Drug name');
for i=1:length(metGenes)
    fprintf(fid,';%s p-value;%s FDR;%s stat;%s identity;%s consumed drug',...
        metGenes{i},metGenes{i},metGenes{i},metGenes{i},metGenes{i});
end
fprintf(fid,';Gene combination p-value;Gene combination FDR;Gene combination stat;Gene combination identity;Gene combination consumed drug');
fprintf(fid,'\n');
for i=1:length(printDrugs)
    fprintf(fid,'%s',drug_names_unique{printDrugs(i)});
    for j=1:size(metGenesDrugs_enrFDR,1)
        if ~isempty(metGenesDrugs_enrStat{j,printDrugs(i)})
            fprintf(fid,';%.8f;%.8f;(%d,%d,%d,%d);%.3f;%.3f',...
                metGenesDrugs_enrP(j,printDrugs(i)),...
                metGenesDrugs_enrFDR(j,printDrugs(i)),...
                metGenesDrugs_enrStat{j,printDrugs(i)},...
                metGenesDrugs_enrPident(j,printDrugs(i)),...
                metGenesDrugs_enrFC(j,printDrugs(i)));
        else
             fprintf(fid,';%.8f;%.8f;(0,0,0,0);%.3f;%.3f',...
                metGenesDrugs_enrP(j,printDrugs(i)),...
                metGenesDrugs_enrFDR(j,printDrugs(i)),...
                metGenesDrugs_enrPident(j,printDrugs(i)),...
                metGenesDrugs_enrFC(j,printDrugs(i)));
        end
    end
    if ~isempty(metGenesDrugs_combined_enrStat{printDrugs(i)})
        fprintf(fid,';%.8f;%.8f;(%d,%d,%d,%d);%.3f;%.3f',...
                metGenesDrugs_combined_enrP(printDrugs(i)),...
                metGenesDrugs_combined_enrFDR(printDrugs(i)),...
                metGenesDrugs_combined_enrStat{printDrugs(i)},...
                metGenesDrugs_combined_enrPident(printDrugs(i)),...
                metGenesDrugs_combined_enrFC(printDrugs(i)));
    else
        fprintf(fid,';%.8f;%.8f;(0,0,0,0);%.3f;%.3f',...
                metGenesDrugs_combined_(printDrugs(i)),...
                metGenesDrugs_combined_(printDrugs(i)),...
                metGenesDrugs_combined_(printDrugs(i)),...
                metGenesDrugs_combined_(printDrugs(i)));
    end
        
    fprintf(fid,'\n');
end
fclose(fid);

