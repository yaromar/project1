
# coding: utf-8

# In[17]:


#!/usr/bin/env python 3
import re

PATH = "./"
#PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "chains.tsv"
PDB_LIST = "pdbList.tsv"


# In[18]:


#PARAMETERS:
#fils is a string with path to the file to be checked

#RESULTS:
#Returns 1 or 0 depending on file existence 


def file_check(file):
    
    try:
        open(file, "r")
        return 1
    
    except IOError:
        #print("Error: " + file + " does not appear to exist.")
        return 0


# In[19]:


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
                typeCount[0] += 'h2a|'
                
            elif(re.search(r'h2b', name, re.I)):
                typeCount[0] += 'h2b|'
                    
            elif(re.search(r'h3', name, re.I)):
                typeCount[0] += 'h3|'
                    
            elif(re.search(r'h4', name, re.I)):
                typeCount[0] += 'h4|'
                    
            elif(re.search(r'h1', name, re.I)):
                typeCount[0] += 'h1|'  
                    
            elif(re.search(r'h5', name, re.I)):
                typeCount[0] += 'h5|'
                    
            else:
                typeCount[0] += 'some histone|'


# In[20]:


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


# In[21]:


#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS: 
#The format of the end-product dictionary is: {pdb : {AlexChain: myChain|UNIPROT|name|type|process|function|organism|nucleosome(bool)|bindingPartner(bool)}}
#Example: {1alq : {'G': 'E|P02302|Histone H3.3C|H3|1212|4141|human|1|0'}}


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
                #print("Error: " + file + " does not appear to exist.")
        
        
        
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
                
                #######################
                if(tempCount): #if the chain is a [part of a] histone
                    
                    if(pdb in histoneCount):

                        if(tempType not in histoneCount[pdb]):
                            histoneCount[pdb].append(tempType) #!!!!!!

                    else:
                        histoneCount[pdb] = [tempType]       
                        
                try: #adds a chain entry to the dict only if there exists a corresponding chain in the mapping file
                    alexChain = list(tempDict[pdb].keys())[list(tempDict[pdb].values()).index(chain)]

                    if(pdb in dictionary):

                        if(tempCount): #checks if chain is a histone!!!!
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType #!!!!

                        else: #!!!!
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|'#!!!!

                    else:

                        if(tempCount): #checks if chain is a histone!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType} #!!!!

                        else: #!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|'} #!!!!
                        
                except ValueError:
                    #print("Error: " + str(ValueError) + ", in " + pdb)               
                    pass

                    
                    
        for structure in list(dictionary): 
            
            if(structure in histoneCount):
                uniqueHistoneNum = len(histoneCount[structure])
                partnerFlag = 0
            
                if(uniqueHistoneNum > 2): #checks if pdb has at least a 3 different histones ~ is a nucleosome!!!!!!!

                    for chain in dictionary[structure]:
                        chainFields = dictionary[structure][chain].split('|')

                        chainType = chainFields[3]##

                        dictionary[structure][chain] += 'nucleosome:1|' #!!!!!!

                        if(partnerFlag == 0 and chainType == 'other'):
                            partnerFlag = 1

                    if(partnerFlag == 0):

                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:0|'

                    else:

                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1|'

                else: #!!!!!!

                    for chain in dictionary[structure]: #!!!!!
                        chainFields = dictionary[structure][chain].split('|')

                        chainType = chainFields[3]

                        dictionary[structure][chain] += 'nucleosome:0|' #!!!!!! 

                        if(partnerFlag == 0 and chainType == 'other'):
                            partnerFlag = 1

                    if(partnerFlag == 0):

                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:0|'
                            
                    else:

                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1|'
            
            else:
                del dictionary[structure]


# In[22]:


def residue_count(interfaceFiles, chainDictionary, interfaceDictionary):
    
    for file in interfaceFiles:
        pdb = file.split('/')[-1].split('_', 1)[0] 
        
        try:
            
            with open (file, 'r') as ifh:
                with open ('hitdata.txt', 'r') as hitFile:
                    ifh.readline()
                
                    for line in ifh:
                        lineFields = line.split('\t')

                        chain1 = lineFields[0].split('_', 1)[0] # the split part treats biological assembly chains as separate chains ???
                        chain2 = lineFields[4].split('_', 1)[0]

                        residue1 = lineFields[2]
                        residue2 = lineFields[6]

                        fields1 = chainDictionary[pdb][chain1].split('|')
                        fields2 = chainDictionary[pdb][chain2].split('|')

                        uniprot1 = fields1[1]
                        uniprot2 = fields2[1]
         
                        type1 = fields1[3]
                        type2 = fields2[3]
                  
                        if(type1 == 'other' and type2 != 'other'):
                            for hitLine in hitFile:
                                li = hitLine.strip()
                                if(li.startswith('Q#')):
                                    lineFields = li.split('\t')
                                    if(lineFields[0].split(' - ')[1] == uniprot1):
                                        if(int(residue1) >= int(lineFields[3]) and int(residue1) <= int(lineFields[4])):
                                            print(type2 + '\t' + lineFields[8])
                            hitFile.seek(0)
                       
                        elif(type2 == 'other' and type1 != 'other'):
                            for hitLine in hitFile:
                                li = hitLine.strip()
                                if(li.startswith('Q#')):
                                    lineFields = li.split('\t')
                                    if(lineFields[0].split(' - ')[1] == uniprot2):
                                        if(int(residue2) >= int(lineFields[3]) and int(residue2) <= int(lineFields[4])):
                                            print(type1 + '\t' + lineFields[8])
                            hitFile.seek(0)                            
      
                           
                  
        except (IOError, KeyError) as e:
            #print("Error: " + file + " does not appear to exist.")
            pass


# In[23]:


def main():
    chainDictionary = {}
    chainDictionary['chain'] = {}
    get_chain_dictionaries(CHAIN_FILE, chainDictionary)
    #for pdb in chainDictionary:
        #for chain in chainDictionary[pdb]:
            #print(pdb + '\t' + chain + '\t' + str(chainDictionary[pdb][chain]))
            
    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')

    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    residue_count(interfaceFiles, chainDictionary, interfaceDictionary)
    
    


# In[24]:


if __name__ == "__main__":
    main()

