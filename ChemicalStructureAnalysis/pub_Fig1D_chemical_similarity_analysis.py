# -*- coding: utf-8 -*-
"""
% script pub_Fig1D_chemical_similarity_analysis
% perform chemical similarity analysis of the defined drug clusters
% with RDkit tools
% plot structures with highlighted common substructures and save to png file

Requiremens: RDKit, pubchempy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
"""

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
import pubchempy as pcp


#############################################################################
filenames = ['Fig1D_chemical_similarity_cluster1.png',
             'Fig1D_chemical_similarity_cluster2a.png',
             'Fig1D_chemical_similarity_cluster2b.png']
#############################################################################
            # Drugs from cluster I from figure 1C
            # (degraded mainly by Bacteroidetes)
clustermolnames = [['Racecadotril','Norethindrone acetate','Roxatidine acetate',
           'Famciclovir','Vilazodone'],
            # Drugs from cluster II from figure 1C
            # (degraded by most bacteria but Proteobacteria)
            ['Tinidazole', 'Entacapone', 'Nitrendipine'],
            ['Phenazopyridine', 'Sulfasalazine']]

# loop through the defined chemical clusters
for i in range(len(filenames)):
    molnames = clustermolnames[i]
    mymols = []
    # loop throug hthe molecules in current cluster
    for molname in molnames:
        #load molecule from PubChem
        results = pcp.get_compounds(molname, 'name')
        mycompound = results[0]
        mymol = Chem.MolFromSmiles(mycompound.canonical_smiles)
        mymols.append(mymol)
    
    # find the maximum common substructure of the structures 
    # forbid matching of ring parts to non-ring parts
    myxcommol = rdFMCS.FindMCS(mymols, ringMatchesRingOnly=True)
    # convert the structure to chem molecule from smarts
    myxcommol = Chem.MolFromSmarts(myxcommol.smartsString)
    # draw a grid highlighting the common substructure
    img = Draw.MolsToGridImage( mymols,
                                subImgSize=(400, 400),
                                legends = molnames,
                                highlightAtomLists=[ mol.GetSubstructMatch(myxcommol) for mol in mymols],
                                useSVG=False )
    # save image to file
    img.save(filenames[i])
