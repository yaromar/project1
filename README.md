# PROJECT1
Characterizing Binding Partners/Spots of Nucleosomes


FILES:

  notes_Y.txt
  
    *tips/useful links/protocols

  nucleosomes_pdb:
 
    *contains PDB# of all (16/7/18) resolved whole (8 chains) nucleosomes and a structure of a chromatosome; 159 structures in total              
    *1 column: PDB
    *'Nucleosome_PDB' entity from the 'ER_diagram'
    *the list can be obtained by completing the first step from 'notes_Y.txt'
    
  protein_ligand_merged.csv:
  
    *contains all chains/nucleotides/ligands of the nucleosome structures from 'nucleosome_pdb'
    *15 columns: structureId, chainId, structureTitle, pdbDoi, entityId, ligandId, ligandName, InChI, InChIKey, hetId, IC50, deltaG, uniprotAcc, uniprotRecommendedName, uniprotAlternativeNames
    *merged 'Histone' and 'Molecule' entities from 'ER_diagram
    *the file can be obtained by completing the second step from 'notes_Y.txt'(additional info can be added)
      
  non_histone_chains.csv:
  
    *contains UNIPROT# of chains that are not histones found in the PDB structures from 'nucleosome_pdb'
    *1 column: UNIPROT
    *there are 30 chains in total, but it has to be double-checked, as it was done manually 
    *some of them might not bind to nucleosome core particles (e.g., bind to DNA)
    
  small_molecules.txt:
  
    *contains ligand names of small molecules found in the PDB structures from 'nucleosome_pdb'
    *1 column: ligandName
    *there are 31 small molecules in total, but it has to be double double-checked, as it was done manually 
    *some of them might not bind to nucleosome core particles (e.g., bind to other chains)
  
  report_1.pdf:
  
    *names of unique binding partner chains extracted from 'protein_ligand_merged.csv' (manually! 30 in total)
    *names of unique small molecules extracted from 'protein_ligand_merged.csv' (manually! omitting ions! 30 in total; 1 missing in the report)
    *tables are from the PIR server, pie charts from the PantherDB 
    *PantherDB recognized only 19 UNIPROT
    
SCRIPTS:

  getpdb.sh:
  
    *downloads a specified PDB to a specified folder (downloads to the current folder by default)
    *usage: getpdb 'PDB#' [directory]
    
