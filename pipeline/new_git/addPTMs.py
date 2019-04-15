
# coding: utf-8

# In[20]:


import csv


# In[24]:


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
            
            for ptmLine in ptmFile:
                if(ptmLine[3] == pdb and ptmLine[1] == uniprot and ptmLine[2] == position):
                    output += '\t' + ptmLine[7] + ' ' + position + ' site'
                    modification = position
                elif(ptmLine[3] == pdb and ptmLine[1] == uniprot):
                    output += '\t' + ptmLine[7] + ' ' + position
                    
            print(output)
            
            ptmfh.seek(0)
        
        else:
            print(output)

