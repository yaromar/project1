# PROJECT1
Characterizing Binding Partners/Spots of Nucleosomes


FILES:

  notes_Y.txt
  
    *tips/useful links/protocols
    
  protein_ligand_merged.csv:
  
    *contains all chains/nucleotides/ligands of the nucleosome structures (175 structures)
    *15 columns: structureId, chainId, structureTitle, pdbDoi, entityId, ligandId, ligandName, InChI, InChIKey, hetId, IC50, 
     deltaG, uniprotAcc, uniprotRecommendedName, uniprotAlternativeNames
    *merged 'Histone' and 'Molecule' entities from 'ER_diagram
    *the file can be obtained by completing the second step from 'notes_Y.txt'(additional info can be added)
    *5OY7 is a tetranucleosome 
    *5T5K is not on the list as it's only 6 chains
    
  non_histone_chains.tsv:
  
    *contains UNIPROT# of chains that are not histones found in the PDB structures from 'nucleosome_pdb'
    *9 columns: UNIPROT,	PDB,	binding_partner_name,	PDB_name,	nucleosome_Organism,	binding_parner_organism, notes, 
     CD_results, deposit_year, release_year
    *there are 31 unique chains in total (if we count centromere protein N in 6C0W and 6BUZ), but 38 entries, as some of the 
     binding partners occur in several PDB structures; also, 6FML's chain G is actually a catalytic domain of INO80 (Ino80), 
     but it does not have a UNIPROT -- I distinguish it from INO80 complex subunit B of 6ETX for now (it aligned with RMSD 
     3.46 2646 atoms)
    *there are only 20 PDB structures that contain binding partners; also, 6C0W is bound only to centromere protein N, 
     6BUZ's binding partner is the same chain as centromere protein N, and 5x0x does not have an interface involving that
     involves histones with binding partners
      
  report_1.pdf:
  
    *names of unique binding partner chains extracted from 'protein_ligand_merged.csv' (manually! 31 in total; some are 
    missing in 
     the report)
    *names of unique small molecules extracted from 'protein_ligand_merged.csv' (manually! omitting ions! some 
     are missing in the report)
    *tables are from the PIR server, pie charts from the PantherDB 
    *PantherDB recognized only 19 UNIPROT
    
  Interfaces/:
  
    *contains results of Alexander's pipeline for PDB structures that contain binding partners
    
  PDB_release_trend*.tsv:
  
    *contains data on years/pdb_structures for histograms in Images/'PDB_release_trends.png'
    *6 columns: year, nucleosome_w/_bp, nucleosome_w/bp_cum, nucleosome_only, nucleosome_only_cum, all_nucleosomes
    
  labeled_chains.tsv:
  
    *classifies chains from "non_histone_chains.tsv" by categories
    *4 columns: PDB, HISTONES_CHAINS, BP_CHAINS, DNA_CHAINS
    
  aa_freqByPDB.tsv:
  
    *contains count of interface residues sorted by pdb and then by chain
    *4 columns: pdb, chain, residue, count
    *does not have 5x0x interface as the closest histone/bp residues are further than 5A away from each other
    
  aa_freqByChain.tsv:
  
    *contains count of interface residues sorted by chains (sum of all PDBs)
    *3 columns: chain, residue, count
    
    
SCRIPTS:

  getpdb.sh:
  
    *downloads a specified PDB to a specified folder (downloads to the current folder by default)
    *usage: getpdb 'PDB#' [directory]
    
  interface_to_frequency.ipynb:
  
    *converts interface data obtained by Alex to residue frequency
    *uses *chain_protein_mapping.tab files and labeled_chains.tsv
    *usage: to be determined

  freqByPDB_to_freqByChain.ipynb
  
    *converts data obtained from interface_to_frequency.ipynb to residue count by chains
    *uses aa_freqByPDB.tsv
    *usage: to be determined
    *does not use 6fml and 6etx as they are cryo em structures, and chains cannot be assigned to canonical nucleosome   
     histones
