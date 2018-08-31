
# coding: utf-8

# In[3]:



#NOTES:
#To get a list of PDB that contains histones from 'text.csv':
#cut -f1 text.tsv | uniq | awk '{print tolower($0)}' | sort

#text.tsv should be sorted by pdb and then uniprot name
#it MUST have 'NA' in uniprot and name blanks


# In[4]:


#!/usr/bin/env python 3
import re

PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "text.tsv"
PDB_LIST = "pdbList.txt"


# In[5]:


#PARAMETERS:
#fils is a string with path to the file to be checked

#RESULTS:
#Returns 1 or 0 depending on file existence 


def file_check(file):
    
    try:
        open(file, "r")
        return 1
    
    except IOError:
        print("Error: " + file + " does not appear to exist.")
        return 0


# In[6]:


#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone(name, typeCount):

    if(not re.search(r'chaperone|ase|binding|p53 peptide|non-histone|jmjc|rna|synth', name, re.I)):
        if(re.search(r'histone.*h?\d|h?\d.*histone|h?\d.*histone-like|histone-like.*h?\d|histone macro.*h?\d|h?\d.*histone macro|h?\d.*\speptide|\speptide.*h?\d|h3k4me0|h3(1-9)k4me3|$h\d^|archaeal histone|histone peptide', name, re.I)):
            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain)

            if(re.search(r'h2a', name, re.I)):
                typeCount[0] += 'h2a#'
                
            elif(re.search(r'h2b', name, re.I)):
                typeCount[0] += 'h2b#'
                    
            elif(re.search(r'h3', name, re.I)):
                typeCount[0] += 'h3#'
                    
            elif(re.search(r'h4', name, re.I)):
                typeCount[0] += 'h4#'
                    
            elif(re.search(r'h1', name, re.I)):
                typeCount[0] += 'h1#'  
                    
            elif(re.search(r'h5', name, re.I)):
                typeCount[0] += 'h5#'
                    
            else:
                typeCount[0] += 'some histone'


# In[7]:


#PARAMETERS: 
#pdbList is a text file with a header and one column PDB
#files is a list
#parameter is a string, either 'mapping' or 'interface' depending on desired results

#RESULTS:
#A list of absolute paths to either mapping files or interface files as it is stored on local NCBI machines


def get_files(pdbList, files, parameter):
    
    with open(pdbList, 'r') as pfh:
        pfh.readline()
        
        if(parameter == 'mapping'):
            
            for line in pfh:
                line = line.strip()
                folder = line[1] + line[2] 
                files.append(PATH + folder + '/' + line + '_chain_protein_mapping.tab')
                
        elif(parameter == 'interface'): 
            
            for line in pfh:
                line = line.strip()
                folder = line[1] + line[2]
                files.append(PATH + folder + '/' + line + '_atomic_contacts_5.0A.tab')


# In[8]:


#PARAMETERS: 
#pdb is a string with pdb id
#parameter is a string, either 'mapping' or 'interface' depending on desired results

#RESULTS:
#An absolute paths to either mapping file or interface file as it is stored on local NCBI machines


def get_file(pdb, parameter):
    
    if(parameter == 'mapping'):
        folder = pdb[1] + pdb[2] 
        file = (PATH + folder + '/' + pdb + '_chain_protein_mapping.tab')
        return file  
        
    elif(parameter == 'interface'): 
        folder = pdb[1] + pdb[2]
        file = (PATH + folder + '/' + pdb + '_atomic_contacts_5.0A.tab')
        return file


# In[18]:


#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS:
#The format of the end-product dictionary is: {pdb : {AlexChain: myChain#UNIPROT#name#type#nucleosome(bool)}}
#Example: {1alq : {'G': 'E#P02302#Histone H3.3C#H3#1'}}


def get_chain_dictionaries(cFile, dictionary): 
    
    with open(cFile, 'r') as cfh:
        cfh.readline()
        
        histoneCount = {} #is used to count number of histones in a structure!!!!!!!
        
        tempDict = {}
        tempDict['pdb'] = {}
            
        mappingFiles = [] 
        get_files(PDB_LIST, mappingFiles, 'mapping')
        
        
        
        
        for file in mappingFiles:
            try: #adds a pdb entry to the dict only if mapping file exists

                with open(file, 'r') as mfh:
                    mfh.readline() #skips header
                    pdb = file.split('/')[-1].split('_', 1)[0]

                    for mLine in mfh:
                        chainPair = mLine.split('\t', 2)
                        alexChain = chainPair[0]
                        myChain = chainPair[1]
                        
                        if(pdb in tempDict):
                            tempDict[pdb][alexChain] = myChain
                        else:
                            tempDict[pdb] = {alexChain : myChain}

            except IOError:
                pass
                #print("Error: " + mappingFile + " does not appear to exist.")
        
        
        
        
        for cLine in cfh:           
            fields = cLine.strip().split('\t')
            
            pdb = fields[0]
    
            if(pdb in tempDict): #continues only if a mapping file exists
                chain = fields[1]
                uniprot = fields[2]
                name = fields[3]


                histoneTypeAndCount = ['', 0]

                is_histone(name, histoneTypeAndCount) #checks whether the name looks like a histone!!!!!!!!!

                tempType = histoneTypeAndCount[0]
                tempCount = histoneTypeAndCount[1]

                if(pdb in histoneCount):
                    histoneCount[pdb] += tempCount #!!!!!!
                else:
                    histoneCount[pdb] = tempCount

                        
                        
                        
                try: #adds a chain entry to the dict only if there exists a corresponding chain in the mapping file
                    alexChain = list(tempDict[pdb].keys())[list(tempDict[pdb].values()).index(chain)]

                    if(pdb in dictionary):

                        if(tempCount): #checks if chain is a histone!!!!
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '#' + uniprot + '#' + name + '#' + tempType #!!!!

                        else: #!!!!
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '#' + uniprot + '#' + name + '#' + 'other#' #!!!!

                    else:

                        if(tempCount): #checks if chain is a histone!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '#' + uniprot + '#' + name + '#' + tempType} #!!!!

                        else: #!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '#' + uniprot + '#' + name + '#' + 'other#'} #!!!!
                        
                except ValueError:
                    #print("Error: " + ValueError + ", in " + pdb)               
                    pass

                    
                    
                    
        for structure in histoneCount:

            if(histoneCount[structure] > 3): #checks if pdb has at least a half of nucleosome!!!!!!!

                for chain in dictionary[structure]: #!!!!!!
                    dictionary[structure][chain] = dictionary[structure][chain] + '1' #!!!!!!

            else: #!!!!!!

                for chain in dictionary[structure]: #!!!!!
                    dictionary[structure][chain] = dictionary[structure][chain] + '0' #!!!!!! 


# In[19]:


#PARAMETERS:
#interfaceFiles is a list of strings containing names of interface files
#chainDictionary is a dictionary produced by 'get_chain_dictionaries'
#interfaceDictonary is a nested dictionary with the innermost dict being interfaceDictionary['uniprotPair'] = {}
#Example: 'P84233@P62799': {'44': ['$A#P84233#Histone H3.2#h3#1@B#P62799#Histone H4#h4#1$E#P84233#Histone H3.2#h3#1@F#P62799#Histone H4#h4#1', 17]
#note that data of a chain is separated by #!!!!!!!!


def residue_frequency(interfaceFiles, chainDictionary, interfaceDictionary):
    for file in interfaceFiles:
        pdb = file.split('/')[-1].split('_', 1)[0] 
        
        try:
            with open (file, 'r') as ifh:
                ifh.readline() #skips header

                for line in ifh:
                    lineFields = line.split('\t', 6) #gets only the first 8 columns !!!

                    chain1 = lineFields[0].split('_', 1)[0] # split treats biologica assembly chains as separate chains ???
                    chain2 = lineFields[4].split('_', 1)[0]

                    residue1 = lineFields[2]
                    residue2 = lineFields[5]

                    uniprot1 = chainDictionary[pdb][chain1].split('#', 2)[1]
                    uniprot2 = chainDictionary[pdb][chain2].split('#', 2)[1]

                    uniprotPair1 = uniprot1 + '@' + uniprot2
                    uniprotPair2 = uniprot2 + '@' + uniprot1

                    chainPair1 = chainDictionary[pdb][chain1] + '@' + chainDictionary[pdb][chain2]
                    chainPair2 = chainDictionary[pdb][chain2] + '@' + chainDictionary[pdb][chain1]

                    if(uniprotPair1 in interfaceDictionary):

                        if(residue1 in interfaceDictionary[uniprotPair1]):
                            interfaceDictionary[uniprotPair1][residue1][1] += 1  

                            if(chainPair1 not in interfaceDictionary[uniprotPair1][residue1][0]):
                                interfaceDictionary[uniprotPair1][residue1][0] += '$' + chainPair1                       

                        else:
                            interfaceDictionary[uniprotPair1][residue1] = [chainPair1, 1]           

                    elif(uniprotPair2 in interfaceDictionary):

                        if(residue2 in interfaceDictionary[uniprotPair2]):
                            interfaceDictionary[uniprotPair2][residue2][1] += 1

                            if(chainPair2 not in interfaceDictionary[uniprotPair2][residue2][0]):
                                interfaceDictionary[uniprotPair2][residue2][0] += '$' + chainPair2                        

                        else:
                            interfaceDictionary[uniprotPair2][residue2] = [chainPair2, 1]

                    else:
                         interfaceDictionary[uniprotPair1] = {residue1 : [chainPair1, 1]}  
        except IOError:
            #print("Error: " + interfaceFile + " does not appear to exist.")
            pass


# In[20]:


def main():
    print('hmmm')
    chainDictionary = {}
    chainDictionary['chain'] = {}
    
    get_chain_dictionaries(CHAIN_FILE, chainDictionary)
    
    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')
    
    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    
    residue_frequency(interfaceFiles, chainDictionary, interfaceDictionary)
    
    print('wth')
    print(interfaceDictionary)


# In[21]:


if __name__ == "__main__":
    main()

