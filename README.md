# PROJECT1
Characterizing Binding Partners/Spots of Nucleosomes


FILES:
  notes_Y.txt
    *tips/useful links/protocols

  nucleosomes_pdb:
    1) contains PDB# of all (16/7/18) resolved whole (8 chains) nucleosomes; 158 structures in total
    2) 1 column: PDB
    3) the list can be obtained by completing the first step from 'notes_Y.txt'
    
  protein_ligand_merged.csv:
    *contains all chains/nucleotides/ligands of the nucleosome structures from 'nucleosome_pdb'
    *15 columns: structureId, chainId, structureTitle, pdbDoi, entityId, ligandId, ligandName, InChI, InChIKey, hetId, IC50, deltaG, uniprotAcc, uniprotRecommendedName, uniprotAlternativeNames
    *the file can be obtained by completing the second step from 'notes_Y.txt'(additional info can be added)
    

SCRIPTS:
  getpdb.sh 
    *downloads specified PDB to specified folder 
    *usage: getpdb 'PDB#' [directory]

