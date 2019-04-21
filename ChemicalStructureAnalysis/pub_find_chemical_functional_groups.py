# -*- coding: utf-8 -*-
"""
% script pub_find_chemical_functional_groups
% create chemical functional group matrix based on structures
% with RDkit tools
% plot structures with highlighted common substructures and save to png file

Requiremens: RDKit, pubchempy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
"""


from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
import pubchempy as pcp
from rdkit.Chem import Fragments

# input drug file in sdf format
sdfFile = 'Pharmakon_selected_drugs.sdf'
# output file with functional groups per drug
outFileFunctionalGroups = 'pharmacon_functional_groups_and_pubchemCID_all_new.csv'

###############################################
# load molecules from SDF file
framePharmacon = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule', 
                            includeFingerprints=True)

# calculate presence of functional groups available in RDKit
Number_of_aliphatic_carboxylic_acids = [Chem.Fragments.fr_Al_COO(x) for x in framePharmacon["Molecule"]]
Number_of_aliphatic_hydroxyl_groups = [Chem.Fragments.fr_Al_OH(x) for x in framePharmacon["Molecule"]]
Number_of_aliphatic_hydroxyl_groups_excluding_tert_OH = [Chem.Fragments.fr_Al_OH_noTert(x) for x in framePharmacon["Molecule"]]
Number_of_N_functional_groups_attached_to_aromatics = [Chem.Fragments.fr_ArN(x) for x in framePharmacon["Molecule"]]
Number_of_Aromatic_carboxylic_acide = [Chem.Fragments.fr_Ar_COO(x) for x in framePharmacon["Molecule"]]
Number_of_aromatic_nitrogens = [Chem.Fragments.fr_Ar_N(x) for x in framePharmacon["Molecule"]]
Number_of_aromatic_amines = [Chem.Fragments.fr_Ar_NH(x) for x in framePharmacon["Molecule"]]
Number_of_aromatic_hydroxyl_groups = [Chem.Fragments.fr_Ar_OH(x) for x in framePharmacon["Molecule"]]
Number_of_carboxylic_acids = [Chem.Fragments.fr_COO(x) for x in framePharmacon["Molecule"]]
Number_of_carboxylic_acids = [Chem.Fragments.fr_COO2(x) for x in framePharmacon["Molecule"]]
Number_of_carbonyl_O = [Chem.Fragments.fr_C_O(x) for x in framePharmacon["Molecule"]]
Number_of_carbonyl_O__excluding_COOH = [Chem.Fragments.fr_C_O_noCOO(x) for x in framePharmacon["Molecule"]]
Number_of_thiocarbonyl = [Chem.Fragments.fr_C_S(x) for x in framePharmacon["Molecule"]]
Number_of_C_OH_CCN_Ctert_alkyl_or_C_OH_CCNcyclic = [Chem.Fragments.fr_HOCCN(x) for x in framePharmacon["Molecule"]]
Number_of_Imines = [Chem.Fragments.fr_Imine(x) for x in framePharmacon["Molecule"]]
Number_of_Tertiary_amines = [Chem.Fragments.fr_NH0(x) for x in framePharmacon["Molecule"]]
Number_of_Secondary_amines = [Chem.Fragments.fr_NH1(x) for x in framePharmacon["Molecule"]]
Number_of_Primary_amines = [Chem.Fragments.fr_NH2(x) for x in framePharmacon["Molecule"]]
Number_of_hydroxylamine_groups = [Chem.Fragments.fr_N_O(x) for x in framePharmacon["Molecule"]]
Number_of_XCCNR_groups = [Chem.Fragments.fr_Ndealkylation1(x) for x in framePharmacon["Molecule"]]
Number_of_tert_alicyclic_amines__no_heteroatoms__not_quinine_like_bridged_N_ = [Chem.Fragments.fr_Ndealkylation2(x) for x in framePharmacon["Molecule"]]
Number_of_H_pyrrole_nitrogens = [Chem.Fragments.fr_Nhpyrrole(x) for x in framePharmacon["Molecule"]]
Number_of_thiol_groups = [Chem.Fragments.fr_SH(x) for x in framePharmacon["Molecule"]]
Number_of_aldehydes = [Chem.Fragments.fr_aldehyde(x) for x in framePharmacon["Molecule"]]
Number_of_alkyl_carbamates__subject_to_hydrolysis_ = [Chem.Fragments.fr_alkyl_carbamate(x) for x in framePharmacon["Molecule"]]
Number_of_alkyl_halides = [Chem.Fragments.fr_alkyl_halide(x) for x in framePharmacon["Molecule"]]
Number_of_allylic_oxidation_sites_excluding_steroid_dienone = [Chem.Fragments.fr_allylic_oxid(x) for x in framePharmacon["Molecule"]]
Number_of_amides = [Chem.Fragments.fr_amide(x) for x in framePharmacon["Molecule"]]
Number_of_amidine_groups = [Chem.Fragments.fr_amidine(x) for x in framePharmacon["Molecule"]]
Number_of_anilines = [Chem.Fragments.fr_aniline(x) for x in framePharmacon["Molecule"]]
Number_of_aryl_methyl_sites_for_hydroxylation = [Chem.Fragments.fr_aryl_methyl(x) for x in framePharmacon["Molecule"]]
Number_of_azide_groups = [Chem.Fragments.fr_azide(x) for x in framePharmacon["Molecule"]]
Number_of_azo_groups = [Chem.Fragments.fr_azo(x) for x in framePharmacon["Molecule"]]
Number_of_barbiturate_groups = [Chem.Fragments.fr_barbitur(x) for x in framePharmacon["Molecule"]]
Number_of_benzene_rings = [Chem.Fragments.fr_benzene(x) for x in framePharmacon["Molecule"]]
Number_of_benzodiazepines_with_no_additional_fused_rings = [Chem.Fragments.fr_benzodiazepine(x) for x in framePharmacon["Molecule"]]
Bicyclic = [Chem.Fragments.fr_bicyclic(x) for x in framePharmacon["Molecule"]]
Number_of_diazo_groups = [Chem.Fragments.fr_diazo(x) for x in framePharmacon["Molecule"]]
Number_of_dihydropyridines = [Chem.Fragments.fr_dihydropyridine(x) for x in framePharmacon["Molecule"]]
Number_of_epoxide_rings = [Chem.Fragments.fr_epoxide(x) for x in framePharmacon["Molecule"]]
Number_of_esters = [Chem.Fragments.fr_ester(x) for x in framePharmacon["Molecule"]]
Number_of_ether_oxygens__including_phenoxy_ = [Chem.Fragments.fr_ether(x) for x in framePharmacon["Molecule"]]
Number_of_furan_rings = [Chem.Fragments.fr_furan(x) for x in framePharmacon["Molecule"]]
Number_of_guanidine_groups = [Chem.Fragments.fr_guanido(x) for x in framePharmacon["Molecule"]]
Number_of_halogens = [Chem.Fragments.fr_halogen(x) for x in framePharmacon["Molecule"]]
Number_of_hydrazine_groups = [Chem.Fragments.fr_hdrzine(x) for x in framePharmacon["Molecule"]]
Number_of_hydrazone_groups = [Chem.Fragments.fr_hdrzone(x) for x in framePharmacon["Molecule"]]
Number_of_imidazole_rings = [Chem.Fragments.fr_imidazole(x) for x in framePharmacon["Molecule"]]
Number_of_imide_groups = [Chem.Fragments.fr_imide(x) for x in framePharmacon["Molecule"]]
Number_of_isocyanates = [Chem.Fragments.fr_isocyan(x) for x in framePharmacon["Molecule"]]
Number_of_isothiocyanates = [Chem.Fragments.fr_isothiocyan(x) for x in framePharmacon["Molecule"]]
Number_of_ketones = [Chem.Fragments.fr_ketone(x) for x in framePharmacon["Molecule"]]
Number_of_ketones_excluding_diaryl__a_b_unsat__dienones__heteroatom_on_Calpha = [Chem.Fragments.fr_ketone_Topliss(x) for x in framePharmacon["Molecule"]]
Number_of_beta_lactams = [Chem.Fragments.fr_lactam(x) for x in framePharmacon["Molecule"]]
Number_of_cyclic_esters__lactones_ = [Chem.Fragments.fr_lactone(x) for x in framePharmacon["Molecule"]]
Number_of_methoxy_groups__OCH3 = [Chem.Fragments.fr_methoxy(x) for x in framePharmacon["Molecule"]]
Number_of_morpholine_rings = [Chem.Fragments.fr_morpholine(x) for x in framePharmacon["Molecule"]]
Number_of_nitriles = [Chem.Fragments.fr_nitrile(x) for x in framePharmacon["Molecule"]]
Number_of_nitro_groups = [Chem.Fragments.fr_nitro(x) for x in framePharmacon["Molecule"]]
Number_of_nitro_benzene_ring_substituents = [Chem.Fragments.fr_nitro_arom(x) for x in framePharmacon["Molecule"]]
Number_of_non_ortho_nitro_benzene_ring_substituents = [Chem.Fragments.fr_nitro_arom_nonortho(x) for x in framePharmacon["Molecule"]]
Number_of_nitroso_groups__excluding_NO2 = [Chem.Fragments.fr_nitroso(x) for x in framePharmacon["Molecule"]]
Number_of_oxazole_rings = [Chem.Fragments.fr_oxazole(x) for x in framePharmacon["Molecule"]]
Number_of_oxime_groups = [Chem.Fragments.fr_oxime(x) for x in framePharmacon["Molecule"]]
Number_of_para_hydroxylation_sites = [Chem.Fragments.fr_para_hydroxylation(x) for x in framePharmacon["Molecule"]]
Number_of_phenols = [Chem.Fragments.fr_phenol(x) for x in framePharmacon["Molecule"]]
Number_of_phenolic_OH_excluding_ortho_intramolecular_Hbond_substituents = [Chem.Fragments.fr_phenol_noOrthoHbond(x) for x in framePharmacon["Molecule"]]
Number_of_phosphoric_acid_groups = [Chem.Fragments.fr_phos_acid(x) for x in framePharmacon["Molecule"]]
Number_of_phosphoric_ester_groups = [Chem.Fragments.fr_phos_ester(x) for x in framePharmacon["Molecule"]]
Number_of_piperdine_rings = [Chem.Fragments.fr_piperdine(x) for x in framePharmacon["Molecule"]]
Number_of_piperzine_rings = [Chem.Fragments.fr_piperzine(x) for x in framePharmacon["Molecule"]]
Number_of_primary_amides = [Chem.Fragments.fr_priamide(x) for x in framePharmacon["Molecule"]]
Number_of_primary_sulfonamides = [Chem.Fragments.fr_prisulfonamd(x) for x in framePharmacon["Molecule"]]
Number_of_pyridine_rings = [Chem.Fragments.fr_pyridine(x) for x in framePharmacon["Molecule"]]
Number_of_quarternary_nitrogens = [Chem.Fragments.fr_quatN(x) for x in framePharmacon["Molecule"]]
Number_of_thioether = [Chem.Fragments.fr_sulfide(x) for x in framePharmacon["Molecule"]]
Number_of_sulfonamides = [Chem.Fragments.fr_sulfonamd(x) for x in framePharmacon["Molecule"]]
Number_of_sulfone_groups = [Chem.Fragments.fr_sulfone(x) for x in framePharmacon["Molecule"]]
Number_of_terminal_acetylenes = [Chem.Fragments.fr_term_acetylene(x) for x in framePharmacon["Molecule"]]
Number_of_tetrazole_rings = [Chem.Fragments.fr_tetrazole(x) for x in framePharmacon["Molecule"]]
Number_of_thiazole_rings = [Chem.Fragments.fr_thiazole(x) for x in framePharmacon["Molecule"]]
Number_of_thiocyanates = [Chem.Fragments.fr_thiocyan(x) for x in framePharmacon["Molecule"]]
Number_of_thiophene_rings = [Chem.Fragments.fr_thiophene(x) for x in framePharmacon["Molecule"]]
Number_of_unbranched_alkanes_of_at_least_4_members__excludes_halogenated_alkanes_ = [Chem.Fragments.fr_unbrch_alkane(x) for x in framePharmacon["Molecule"]]
Number_of_urea_groups = [Chem.Fragments.fr_urea(x) for x in framePharmacon["Molecule"]]

# combine dataframe
combinedFG = [Number_of_aliphatic_carboxylic_acids, Number_of_aliphatic_hydroxyl_groups, Number_of_aliphatic_hydroxyl_groups_excluding_tert_OH, Number_of_N_functional_groups_attached_to_aromatics, Number_of_Aromatic_carboxylic_acide, Number_of_aromatic_nitrogens, Number_of_aromatic_amines, Number_of_aromatic_hydroxyl_groups, Number_of_carboxylic_acids, Number_of_carboxylic_acids, Number_of_carbonyl_O, Number_of_carbonyl_O__excluding_COOH, Number_of_thiocarbonyl, Number_of_C_OH_CCN_Ctert_alkyl_or_C_OH_CCNcyclic, Number_of_Imines, Number_of_Tertiary_amines, Number_of_Secondary_amines, Number_of_Primary_amines, Number_of_hydroxylamine_groups, Number_of_XCCNR_groups, Number_of_tert_alicyclic_amines__no_heteroatoms__not_quinine_like_bridged_N_, Number_of_H_pyrrole_nitrogens, Number_of_thiol_groups, Number_of_aldehydes, Number_of_alkyl_carbamates__subject_to_hydrolysis_, Number_of_alkyl_halides, Number_of_allylic_oxidation_sites_excluding_steroid_dienone, Number_of_amides, Number_of_amidine_groups, Number_of_anilines, Number_of_aryl_methyl_sites_for_hydroxylation, Number_of_azide_groups, Number_of_azo_groups, Number_of_barbiturate_groups, Number_of_benzene_rings, Number_of_benzodiazepines_with_no_additional_fused_rings, Bicyclic, Number_of_diazo_groups, Number_of_dihydropyridines, Number_of_epoxide_rings, Number_of_esters, Number_of_ether_oxygens__including_phenoxy_, Number_of_furan_rings, Number_of_guanidine_groups, Number_of_halogens, Number_of_hydrazine_groups, Number_of_hydrazone_groups, Number_of_imidazole_rings, Number_of_imide_groups, Number_of_isocyanates, Number_of_isothiocyanates, Number_of_ketones, Number_of_ketones_excluding_diaryl__a_b_unsat__dienones__heteroatom_on_Calpha, Number_of_beta_lactams, Number_of_cyclic_esters__lactones_, Number_of_methoxy_groups__OCH3, Number_of_morpholine_rings, Number_of_nitriles, Number_of_nitro_groups, Number_of_nitro_benzene_ring_substituents, Number_of_non_ortho_nitro_benzene_ring_substituents, Number_of_nitroso_groups__excluding_NO2, Number_of_oxazole_rings, Number_of_oxime_groups, Number_of_para_hydroxylation_sites, Number_of_phenols, Number_of_phenolic_OH_excluding_ortho_intramolecular_Hbond_substituents, Number_of_phosphoric_acid_groups, Number_of_phosphoric_ester_groups, Number_of_piperdine_rings, Number_of_piperzine_rings, Number_of_primary_amides, Number_of_primary_sulfonamides, Number_of_pyridine_rings, Number_of_quarternary_nitrogens, Number_of_thioether, Number_of_sulfonamides, Number_of_sulfone_groups, Number_of_terminal_acetylenes, Number_of_tetrazole_rings, Number_of_thiazole_rings, Number_of_thiocyanates, Number_of_thiophene_rings, Number_of_unbranched_alkanes_of_at_least_4_members__excludes_halogenated_alkanes_, Number_of_urea_groups];
combinedFG = list(map(list, zip(*combinedFG)))

fgDataframe = pd.DataFrame(combinedFG, columns = ["Number_of_aliphatic_carboxylic_acids","Number_of_aliphatic_hydroxyl_groups","Number_of_aliphatic_hydroxyl_groups_excluding_tert_OH","Number_of_N_functional_groups_attached_to_aromatics","Number_of_Aromatic_carboxylic_acide","Number_of_aromatic_nitrogens","Number_of_aromatic_amines","Number_of_aromatic_hydroxyl_groups","Number_of_carboxylic_acids","Number_of_carboxylic_acids","Number_of_carbonyl_O","Number_of_carbonyl_O__excluding_COOH","Number_of_thiocarbonyl","Number_of_C_OH_CCN_Ctert_alkyl_or_C_OH_CCNcyclic","Number_of_Imines","Number_of_Tertiary_amines","Number_of_Secondary_amines","Number_of_Primary_amines","Number_of_hydroxylamine_groups","Number_of_XCCNR_groups","Number_of_tert_alicyclic_amines__no_heteroatoms__not_quinine_like_bridged_N_","Number_of_H_pyrrole_nitrogens","Number_of_thiol_groups","Number_of_aldehydes","Number_of_alkyl_carbamates__subject_to_hydrolysis_","Number_of_alkyl_halides","Number_of_allylic_oxidation_sites_excluding_steroid_dienone","Number_of_amides","Number_of_amidine_groups","Number_of_anilines","Number_of_aryl_methyl_sites_for_hydroxylation","Number_of_azide_groups","Number_of_azo_groups","Number_of_barbiturate_groups","Number_of_benzene_rings","Number_of_benzodiazepines_with_no_additional_fused_rings","Bicyclic","Number_of_diazo_groups","Number_of_dihydropyridines","Number_of_epoxide_rings","Number_of_esters","Number_of_ether_oxygens__including_phenoxy_","Number_of_furan_rings","Number_of_guanidine_groups","Number_of_halogens","Number_of_hydrazine_groups","Number_of_hydrazone_groups","Number_of_imidazole_rings","Number_of_imide_groups","Number_of_isocyanates","Number_of_isothiocyanates","Number_of_ketones","Number_of_ketones_excluding_diaryl__a_b_unsat__dienones__heteroatom_on_Calpha","Number_of_beta_lactams","Number_of_cyclic_esters__lactones_","Number_of_methoxy_groups__OCH3","Number_of_morpholine_rings","Number_of_nitriles","Number_of_nitro_groups","Number_of_nitro_benzene_ring_substituents","Number_of_non_ortho_nitro_benzene_ring_substituents","Number_of_nitroso_groups__excluding_NO2","Number_of_oxazole_rings","Number_of_oxime_groups","Number_of_para_hydroxylation_sites","Number_of_phenols","Number_of_phenolic_OH_excluding_ortho_intramolecular_Hbond_substituents","Number_of_phosphoric_acid_groups","Number_of_phosphoric_ester_groups","Number_of_piperdine_rings","Number_of_piperzine_rings","Number_of_primary_amides","Number_of_primary_sulfonamides","Number_of_pyridine_rings","Number_of_quarternary_nitrogens","Number_of_thioether","Number_of_sulfonamides","Number_of_sulfone_groups","Number_of_terminal_acetylenes","Number_of_tetrazole_rings","Number_of_thiazole_rings","Number_of_thiocyanates","Number_of_thiophene_rings","Number_of_unbranched_alkanes_of_at_least_4_members__excludes_halogenated_alkanes_","Number_of_urea_groups"])
framePharmacon_combined = framePharmacon.copy()
framePharmacon_combined = framePharmacon_combined.drop(['Molecule'], axis=1)
framePharmacon_combined = framePharmacon_combined.reset_index()
framePharmacon_combined = framePharmacon_combined.join(fgDataframe)

# find pubchemid and exact mass
pubchem_cid = []
exact_masses = []
for molname in framePharmacon['MOLENAME']:
    results = pcp.get_compounds(molname, 'name')
    if results:
        mycompound = results[0]
        pubchem_cid.append(mycompound.cid)
        exact_masses.append(mycompound.exact_mass)
    else:
        pubchem_cid.append('')
        exact_masses.append(0)
   
fgDataframePubchem = {"PubchemCID":pubchem_cid, "Exact mass":exact_masses}
fgDataframePubchem = pd.DataFrame(fgDataframePubchem)
framePharmacon_combined = framePharmacon_combined.join(fgDataframePubchem)

# save new dataframe with added columns to file
framePharmacon_combined.to_csv(outFileFunctionalGroups)
