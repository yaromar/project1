
# coding: utf-8

# In[3]:


import csv


# In[10]:


with open("combined/PTMs.tsv",'r') as ptmfh, open("combined/uniqueResults.tsv", 'r') as interfh, open("combined/interactome.tsv", 'w') as resfh:
    ptmFile = csv.reader(ptmfh, delimiter =';')
    interFile = csv.reader(interfh, delimiter = '\t')
    
    modDict = {}
    
    for interLine in interFile:
        output = '\t'.join(interLine)
        
        if(interLine[-3] != "other"):
            pdb = interLine[-1].upper()
            uniprot = interLine[0].upper()
            position = interLine[-2]
            
            flagSite = 0
            flagMod = 0
            modification = ''
            for ptmLine in ptmFile:
                if(ptmLine[3] == pdb and ptmLine[1] == uniprot and ptmLine[2] == position and not flagSite and not flagMod):
                    output += '\t' + ptmLine[7] + ' ' + position
                    modification = ptmLine[7]
                    flagSite = 1
                    flagMod = 1
                elif(ptmLine[3] == pdb and ptmLine[1] == uniprot and ptmLine[2] == position and not flagSite):
                    output += '|' + ptmLine[7] + ' ' + position
                    modification = ptmLine[7]
                    flagSite = 1
                    flagMod = 1        
                elif(ptmLine[3] == pdb and ptmLine[1] == uniprot and flagMod):
                    output += '|' + ptmLine[7] + ' ' + ptmLine[2]
                    flagMod = 1
                elif(ptmLine[3] == pdb and ptmLine[1] == uniprot):
                    output += '\t' + ptmLine[7] + ' ' + ptmLine[2]
                    flagMod = 1
                    
            if(flagSite):
                output += '\t' + modification + '\n'
            else:
                output += '\t' + 'NA' + '\t' + 'NA' + '\n'
                
            ptmfh.seek(0)
        
            resfh.write(output)
        else:
            output += '\t' + 'NA' + '\t' + 'NA' + '\n'
            resfh.write(output)

