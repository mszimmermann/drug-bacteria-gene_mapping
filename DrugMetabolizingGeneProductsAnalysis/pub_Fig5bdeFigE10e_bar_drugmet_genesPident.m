%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig5bdeFigE10e_bar_drugmet_genesPident
% read drug metabolism data and gene-genomes identity BLAST file
% plot bar plots of species metabolizing the drug 
% and corresponding identified metabolism gene product gene identity
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
% pidentThreshold = 50; % threshold of protein identity for plotting
% outfile_Fig5deFigE10e_bar_drugmetabolism_genes_pident

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
% reformat strain names
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read individual data for drug fold changes
drugFoldChanges_individual = readtable(outfile_Table_individual_drug_fold_changes,...
                            'HeaderLines', 1);
drugAnalysisColumnNames_ind = drugFoldChanges_individual.Properties.VariableNames;
FCcolumns_ind = cellfun(@(x) contains(x, 'x_Consumed'),...
                         drugAnalysisColumnNames_ind);
drugFoldChanges_individual_data = drugFoldChanges_individual{:,FCcolumns_ind};
drugFC12to0_Species_ind = drugAnalysisColumnNames_ind(FCcolumns_ind);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedA', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedB', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedC', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedD', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = reshape(drugFC12to0_Species_ind,4,[])';

drugFC12to0_Species_ind = cellfun(@(x) strrep(x, '_', ''), drugFC12to0_Species_ind, 'unif', 0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x, '-', ''), drugFC12to0_Species_ind, 'unif', 0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x, '.', ''), drugFC12to0_Species_ind, 'unif', 0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x, '(', ''), drugFC12to0_Species_ind, 'unif', 0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x, ')', ''), drugFC12to0_Species_ind, 'unif', 0);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'reuteri'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'reuteri'),AssaySpecies(:,4)), 1});
drugFC12to0_SpeciesS{cellfun(@(x) contains(x, 'WH2'),drugFC12to0_Species)} = ...
    sprintf('S%03d', AssaySpecies{cellfun(@(x) contains(x, 'wh2'),AssaySpecies(:,4)), 1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match individual drug data
drugFC12to0_Species_for_ind = drugFC12to0_Species;
drugFC12to0_Species_for_ind = cellfun(@(x) strrep(x, ' ', ''), drugFC12to0_Species_for_ind, 'unif', 0);
[~,~,drug_ind_idx] = intersect(lower(drugFC12to0_Species_for_ind), lower(drugFC12to0_Species_ind(:,1)), 'stable');

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
% for each drug plot species nmetabolism bars
% and species identity to identified drug metabolizing gene products
drug_names_unique = gene_drug_list(:,2);
% remove underscores if there are any in drug names
for i=1:length(drug_names_unique)
    if contains(drug_names_unique{i}, '_')
        drug_names_unique{i} = drug_names_unique{i}(1:strfind(drug_names_unique{i},'_')-1);
    end
end   
drug_names_unique = unique(drug_names_unique);
for drug_i = 1:length(drug_names_unique)
    %find genes that metabolize the drug
    curgeneidx = ismember(gene_drug_list(:,2),drug_names_unique(drug_i));
    
    % combined pident for all drug genes
    drug_clustidx = ismember(lower(drugClustergram.DrugName), lower(drug_names_unique(drug_i)));
    curAdjustedFC = drugFC12to0_t0toCTRL_combined(:, drug_clustidx);
    curAdjustedFCSTD = drugFC12to0_t0toCTRL_STDcombined(:, drug_clustidx);
    
    % get individual drug fold changes
    drugFoldChanges_individual_drugIDX = cellfun(@(x) contains(lower(x), lower(drug_names_unique(drug_i))),...
                                                drugFoldChanges_individual.DrugName);
    cur_individual_data = drugFoldChanges_individual_data(drugFoldChanges_individual_drugIDX,:);
    cur_individual_data = reshape(cur_individual_data,4,[])';
    cur_individual_data = cur_individual_data(drug_ind_idx,:);
    
    if nnz(curgeneidx)
        curgenes = unique(gene_drug_list(curgeneidx,1));
        % remove collinsella genes for famciclovir
        if contains(lower(drug_names_unique(drug_i)), 'famciclovir')
            curgenes(cellfun(@(x) contains(lower(x), 'colaer_'), curgenes)) =[];
        end
        combinedGenePident = zeros(length(curAdjustedFC),length(curgenes));  
         
        
        fig = figure;%('units','normalized','outerposition',[0 0 1 1]);
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
            end
        end
        [~, colidx] = sort(curAdjustedFC);
        % sort fold changes and fold change std according to fold change
        % for plotting
        plotMatrix = combinedGenePident>pidentThreshold;
        
        plotSTD = curAdjustedFCSTD.*(2.^curAdjustedFC)*log(2);
        plotFC = 1-2.^curAdjustedFC;
        
        plotMatrix = plotMatrix(colidx,:);
        plotFC = plotFC(colidx,1);
        plotSTD = plotSTD(colidx,1);
        plot_individual_data = cur_individual_data(colidx,:);
        
        plotMatrix(isnan(combinedGenePident(colidx,1)),:) = [];
        plotFC(isnan(combinedGenePident(colidx,1))) = [];
        plotSTD(isnan(combinedGenePident(colidx,1))) = [];
        plot_individual_data(isnan(combinedGenePident(colidx,1)),:) = [];
        plot_individual_data = plot_individual_data/100;
        plot_individual_data(plot_individual_data==0) = nan;
        % prepare matrix with different color for different gene identity
        plotMatrixGenes = zeros(1,length(curgenes));
        plotMatrixGenes(cellfun(@(x) contains(x, 'BT_'), curgenes)) = 1;
        plotMatrixGenes(cellfun(@(x) contains(x, 'BACDOR_'), curgenes)) = 2;
        plotMatrixGenes(cellfun(@(x) contains(x, 'COLAER_'), curgenes)) = 3;
        plotMatrix = plotMatrix.*repmat(plotMatrixGenes, size(plotMatrix,1),1);
      
        [~, idx] = sort(plotMatrixGenes);
        plotMatrix = plotMatrix(:,idx);
        plotMatrixGenesNames = curgenes(idx);
        % define colormap with three colors for the genes from three species
        mycolormap = [0 0 0;...
                      222 45 38;...
                      255 128 0;...
                      70 170 150;...
                      ]/255;
        genecolors = [1 2 3];
        mycolormap(1+genecolors(~ismember(genecolors, unique(plotMatrixGenes))),:)=[];
        
        % plot bars species drug metabolism 
        subplot(2,1,1)
        bar(1:length(plotFC), plotFC)
        hold on
        errorbar(1:length(plotFC), plotFC, plotSTD, 'k.')
        ylim([0 1.2])
        xlim([0.5 length(plotFC)+0.5])
        hold on
        set(gca, 'YTick', 0:0.2:1)
        set(gca, 'YTickLabel', 0:20:100)
        ylabel('Percent of drug metabolized')
        title(drug_names_unique{drug_i});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add individual data points
        hold on
        scatter(repmat(1:size(plot_individual_data,1),1,4),...
            plot_individual_data(:),...
            20,'b', 'filled')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot species gene pident
        subplot(2,1,2)
        imagesc(plotMatrix')
        hold on
        for i=0.5:1:size(plotMatrix,1)
            plot([i i], [0 length(curgenes)+1], 'k')
        end    
        set(gca, 'YTick', 1:length(plotMatrixGenesNames));
        set(gca, 'YTickLabel', plotMatrixGenesNames);

        set(gca,'TickLabelInterpreter', 'none');
        
        colormap(mycolormap)    
        title(sprintf('Protein pident >=%d', pidentThreshold));
        % save to file
        orient landscape
        print(fig, '-painters', '-dpsc', '-r600', '-bestfit', '-append',...
              outfile_Fig5deFigE10e_bar_drugmetabolism_genes_pident)
        close(fig)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%