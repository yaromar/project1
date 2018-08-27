axes location off

# graphics top delete all

color Display Background white

#
# 1st molecule (also the point of reference)
#

mol load pdb pdb/5x0y.pdb

# Delete default representation
mol delrep 0 top

# mol load pdb pdb/5o9g.pdb

#
# Choose DNA of first molecule for alignment later
#

set sel0 [atomselect top "nucleic and backbone and name P"]
# set sel0 [atomselect top "nucleic and backbone and name P or (alpha and ((chain A and resid 43 to 130) or (chain E and resid 43 to 130)))"]


#
# Visualize
#


#  H3 H3 Dimer representation
mol representation QuickSurf 0.8 0.8 0.5
mol color ColorID 0
mol selection "(chain A or chain E)"
mol material Transparent
mol addrep top
# mol selupdate 0 top 0
# mol colupdate 0 top 0
# mol smoothrep top 0 5


#  H4 H4 Dimer representation
mol representation QuickSurf 0.8 0.8 0.5
mol color ColorID 7
mol selection "chain B or chain F"
mol material Transparent
mol addrep top

#  H2A H2A Dimer representation
mol representation QuickSurf 0.8 0.8 0.5
mol color ColorID 4
mol selection "chain C or chain G"
mol material Transparent
mol addrep top

#  H2B H2B Dimer representation
mol representation QuickSurf 0.8 0.8 0.5
mol color ColorID 1
mol selection "chain D or chain H"
mol material Transparent
mol addrep top

# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 50.000000 2.070000 1
mol color ColorID 6
mol selection "nucleic and backbone"
mol material AOEdgy
mol addrep top

# DNA phosphates in contact
mol representation VDW 1.0 12
mol color ColorID 3
mol selection "nucleic and backbone and name P and within 5 of chain O"
mol material AOShiny
mol addrep top

# Protein atoms in contact
mol representation VDW 1.0 12
mol color ResType
mol selection "(chain A or chain B or chain C or chain D or chain E or chain F or chain G or chain H) and within 5 of chain O"
mol material AOShiny
mol addrep top

#
# 2nd molecule (just visualize interface)
#

mol load pdb pdb/5x0x.pdb

# Delete default representation
mol delrep 0 top

#
# Align wrt 1st molecule
#

set sel1 [atomselect top "nucleic and backbone and name P"]
# set sel1 [atomselect top "alpha and ((chain A and resid 43 to 130) or (chain E and resid 43 to 130))"]
# set sel1 [atomselect top "nucleic and backbone and name P or (alpha and ((chain A and resid 43 to 130) or (chain E and resid 43 to 130)))"]
set sel1all [atomselect top "all"]

set M [measure fit $sel0 $sel1]
$sel1all move $M

# DNA phosphates in contact
mol representation VDW 1.0 12
mol color ColorID 3
mol selection "nucleic and backbone and name P and within 5 of chain O"
mol material AOShiny
mol addrep top

# Protein atoms in contact
mol representation VDW 1.0 12
mol color ResType
mol selection "(chain A or chain B or chain C or chain D or chain E or chain F or chain G or chain H) and within 5 of chain O"
mol material AOShiny
mol addrep 0




