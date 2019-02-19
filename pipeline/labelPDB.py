
# coding: utf-8

# In[38]:


#!/usr/bin/env python 3
import re

CHAIN_FILE = "chains.tsv"


# In[39]:


#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone(name, typeCount):

    if(not re.search(r'chaperone|ase|binding|p53 peptide|non-histone|jmjc|rna|synth', name, re.I)):
        if(re.search(r'histone.*h?\d|h?\d.*histone|h?\d.*histone-like|histone-like.*h?\d|histone macro.*h?\d|h?\d.*histone macro|h?\d.*\speptide|\speptide.*h?\d|h3k4me0|h3(1-9)k4me3|$h\d^|archaeal histone|histone peptide', name, re.I)):
            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain!!!)
            if(re.search(r'h2a', name, re.I)):
                typeCount[0] += 'h2a'
                
            elif(re.search(r'h2b', name, re.I)):
                typeCount[0] += 'h2b'
                    
            elif(re.search(r'h3', name, re.I)):
                typeCount[0] += 'h3'
                    
            elif(re.search(r'h4', name, re.I)):
                typeCount[0] += 'h4'
                    
            elif(re.search(r'h1', name, re.I)):
                typeCount[0] += 'h1'  
                    
            elif(re.search(r'h5', name, re.I)):
                typeCount[0] += 'h5'
                    
            else:
                typeCount[0] += 'some histone'


# In[52]:


def main():
    with open("chains.tsv", 'r') as fh:
        fh.readline()

        histoneCount = {}
        hasBP = {}

        for line in fh:

            fields = line.split('\t')
            pdb = fields[0]
            uniprot = fields[4]
            name = fields[5].strip()
 
            if(uniprot != 'NA'):
                histoneTypeAndCount = ['', 0]

                is_histone(name, histoneTypeAndCount) #checks whether the name looks like a histone

                tempType = histoneTypeAndCount[0]
                tempCount = histoneTypeAndCount[1]

                    #######################
                if(tempCount): #if the chain is a [part of a] histone
                    if(pdb in histoneCount):

                        if(tempType not in histoneCount[pdb]):
                            histoneCount[pdb] += '|'
                            histoneCount[pdb] += tempType #!!!!!!

                    else:
                        histoneCount[pdb] = tempType
                        
                else:

                    if(pdb in hasBP):
                        hasBP[pdb] = 1

                    else:
                        hasBP[pdb] = 1


        for structure in histoneCount:
            if(len(histoneCount[structure].split('|')) > 2):
                if(structure in hasBP):
                    print(structure + '\t' + 'nucleosome' + '\t' + 'has BP')
                else:
                    print(structure + '\t' + 'nucleosome' + '\t' + 'no BP')

            elif(len(histoneCount[structure].split('|')) > 0):
                if(structure in hasBP):
                    print(structure + '\t' + 'histone' + '\t' + 'has BP')
                else:
                    print(structure + '\t' + 'histone' + '\t' + 'no BP')


# In[53]:


if __name__ == "__main__":
    main()

