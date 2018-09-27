
# coding: utf-8

# In[33]:


#NOTES:
#To get a list of PDB that contains histones from 'text.csv':
#cut -f1 text.tsv | uniq | awk '{print tolower($0)}' | sort

#text.tsv should be sorted by pdb and then uniprot name
#it MUST have 'NA' in uniprot and name blanks


# In[34]:


#!/usr/bin/env python 3
import re

#PATH = "../data/Interfaces/"
PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "text.tsv"
PDB_LIST = "pdbList.txt"


# In[35]:


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


# In[36]:


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


# In[37]:


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


# In[38]:


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


# In[39]:


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
                #print("Error: " + mappingFile + " does not appear to exist.")
        
        
        
        for cLine in cfh:           
            fields = cLine.strip().split('\t')

            pdb = fields[0]
            
            if(pdb in tempDict): #continues only if a mapping file exists
                chain = fields[1]
                organism = fields[2]
                process = fields[3]
                function = fields[4]
                uniprot = fields[5]
                name = fields[6]

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
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType + '|' + process + '|' + function + '|' + organism #!!!!

                        else: #!!!!
                            dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|' + '|' + process + '|' + function + '|' + organism#!!!!

                    else:

                        if(tempCount): #checks if chain is a histone!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType + '|' + process + '|' + function + '|' + organism} #!!!!

                        else: #!!!!
                            dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|' + '|' + process + '|' + function + '|' + organism} #!!!!
                        
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


# In[40]:


#PARAMETERS:
#interfaceFiles is a list of strings containing names of interface files
#chainDictionary is a dictionary produced by 'get_chain_dictionaries'
#interfaceDictonary is a nested dictionary with the innermost dict being interfaceDictionary['uniprotPair'] = {}, where the values are lists of the next form [chain pair data, residue number, pdb IDs] 
#Example: 'P84233@P62799': {'44': ['A|P84233|Histone H3.2|h3|1@B|P62799|Histone H4|h4|1$E|P84233|Histone H3.2|h3|1@F|P62799|Histone H4|h4|1', 17, '1zla$q5cl']
#note that data fields of a chain is separated by |, @ separates chain from binding partner chain, $ separates chain pairs and pdb structures


def residue_count(interfaceFiles, chainDictionary, interfaceDictionary):
    
    for file in interfaceFiles:
        pdb = file.split('/')[-1].split('_', 1)[0] 
        
        try:
            
            with open (file, 'r') as ifh:
                ifh.readline()
                
                atomList1 = [] #atomic count
                atomList2 = [] #atomic count   
                
                for line in ifh:
                    lineFields = line.split('\t')

                    chain1 = lineFields[0].split('_', 1)[0] # the split part treats biological assembly chains as separate chains ???
                    chain2 = lineFields[4].split('_', 1)[0]

                    residue1 = lineFields[2]
                    residue2 = lineFields[6]
                    
                    atom1 = lineFields[3] #atomic count
                    atom2 = lineFields[7] #atomic count
                    
                    atomEntry1 = chain1 + residue1 + atom1 #atomic count
                    atomEntry2 = chain2 + residue2 + atom2 #atomic count
                    
                    fields1 = chainDictionary[pdb][chain1].split('|', 5)
                    fields2 = chainDictionary[pdb][chain2].split('|', 5)
                    
                    uniprot1 = fields1[1]
                    uniprot2 = fields2[1]
                    
                    type1 = fields1[3]
                    type2 = fields2[3]
                    
                    uniprotPair1 = uniprot1 + '@' + uniprot2
                    chainPair1 = chainDictionary[pdb][chain1] + '@' + chainDictionary[pdb][chain2]
                    
                    uniprotPair2 = uniprot2 + '@' + uniprot1
                    chainPair2 = chainDictionary[pdb][chain2] + '@' + chainDictionary[pdb][chain1]
                    
                    if(type1 != 'other'):                    
                        
                        if(uniprotPair1 in interfaceDictionary):
                            
                            if(residue1 in interfaceDictionary[uniprotPair1]):
                                
                                if(atomEntry1 not in atomList1): #atomic count
                                    interfaceDictionary[uniprotPair1][residue1][1] += 1  
                                    atomList1.append(atomEntry1) #atomic count
                                    
                                if(chainPair1 not in interfaceDictionary[uniprotPair1][residue1][0]):
                                    interfaceDictionary[uniprotPair1][residue1][0] += '$' + chainPair1                       

                                    if(pdb not in interfaceDictionary[uniprotPair1][residue1][2]):
                                        interfaceDictionary[uniprotPair1][residue1][2] += '$' + pdb

                            else:
                                interfaceDictionary[uniprotPair1][residue1] = [chainPair1, 1, pdb]    #format of the innermost dict      
                                atomList1.append(atomEntry1) #atomic count
                                
                        elif(uniprotPair2 not in interfaceDictionary):
                            interfaceDictionary[uniprotPair1] = {residue1 : [chainPair1, 1, pdb]} 
                            atomList1.append(atomEntry1)#atomic count
                            
                        else:
                                
                            if(residue2 in interfaceDictionary[uniprotPair2]):
                                
                                if(atomEntry2 not in atomList2): #atomic count
                                    interfaceDictionary[uniprotPair2][residue2][1] += 1  
                                    atomList2.append(atomEntry2) #atomic count
                                    
                                if(chainPair2 not in interfaceDictionary[uniprotPair2][residue2][0]):
                                    interfaceDictionary[uniprotPair2][residue2][0] += '$' + chainPair2                      

                                    if(pdb not in interfaceDictionary[uniprotPair2][residue2][2]):
                                        interfaceDictionary[uniprotPair2][residue2][2] += '$' + pdb

                            else:
                                interfaceDictionary[uniprotPair2][residue2] = [chainPair2, 1, pdb]    #format of the innermost dict             
                                atomList2.append(atomEntry2) #atomic count
                                
                    elif(uniprotPair1 not in interfaceDictionary):
                        
                        if(uniprotPair2 in interfaceDictionary):
                            
                            if(residue2 in interfaceDictionary[uniprotPair2]):
                                
                                if(atomEntry2 not in atomList2): #atomic count
                                    interfaceDictionary[uniprotPair2][residue2][1] += 1  
                                    atomList2.append(atomEntry2) #atomic count
                                    
                                if(chainPair2 not in interfaceDictionary[uniprotPair2][residue2][0]):
                                    interfaceDictionary[uniprotPair2][residue2][0] += '$' + chainPair2                       

                                    if(pdb not in interfaceDictionary[uniprotPair2][residue2][2]):
                                        interfaceDictionary[uniprotPair2][residue2][2] += '$' + pdb

                            else:
                                interfaceDictionary[uniprotPair2][residue2] = [chainPair2, 1, pdb]    #format of the innermost dict     
                                atomList2.append(atomEntry2) #atomic count
                                
                        else:
                            interfaceDictionary[uniprotPair2] = {residue2 : [chainPair2, 1, pdb]} 
                            atomList2.append(atomEntry2) #atomic count    
                            
                    else:
                            
                        if(residue1 in interfaceDictionary[uniprotPair1]):
                            
                            if(atomEntry1 not in atomList1): #atomic count
                                interfaceDictionary[uniprotPair1][residue1][1] += 1  
                                atomList1.append(atomEntry1) #atomic count
                                
                            if(chainPair1 not in interfaceDictionary[uniprotPair1][residue1][0]):
                                interfaceDictionary[uniprotPair1][residue1][0] += '$' + chainPair1                       

                                if(pdb not in interfaceDictionary[uniprotPair1][residue1][2]):
                                    interfaceDictionary[uniprotPair1][residue1][2] += '$' + pdb

                        else:
                            interfaceDictionary[uniprotPair1][residue1] = [chainPair1, 1, pdb]    #format of the innermost dict    
                            atomList1.append(atomEntry1) #atomic count
#CHECK FOR REDUNDANCY!!!

                        
#                     if(uniprotPair1 in interfaceDictionary):

#                         if(residue1 in interfaceDictionary[uniprotPair1]):
#                             interfaceDictionary[uniprotPair1][residue1][1] += 1  

#                             if(chainPair1 not in interfaceDictionary[uniprotPair1][residue1][0]):
#                                 interfaceDictionary[uniprotPair1][residue1][0] += '$' + chainPair1                       
                                
#                                 if(pdb not in interfaceDictionary[uniprotPair1][residue1][2]):
#                                     interfaceDictionary[uniprotPair1][residue1][2] += '$' + pdb
                                   
#                         else:
#                             interfaceDictionary[uniprotPair1][residue1] = [chainPair1, 1, pdb]    #format of the innermost dict        

#                     elif(uniprotPair2 in interfaceDictionary):

#                         if(residue2 in interfaceDictionary[uniprotPair2]):
#                             interfaceDictionary[uniprotPair2][residue2][1] += 1

#                             if(chainPair2 not in interfaceDictionary[uniprotPair2][residue2][0]):
#                                 interfaceDictionary[uniprotPair2][residue2][0] += '$' + chainPair2                        

#                                 if(pdb not in interfaceDictionary[uniprotPair2][residue2][2]):
#                                     interfaceDictionary[uniprotPair2][residue2][2] += '$' + pdb
                                    
#                         else:
#                             interfaceDictionary[uniprotPair2][residue2] = [chainPair2, 1, pdb]

#                     else:
#                          interfaceDictionary[uniprotPair1] = {residue1 : [chainPair1, 1, pdb]}  
                            
        except (IOError, KeyError) as e:
            #print("Error: " + interfaceFile + " does not appear to exist.")
            pass


# In[41]:


def normalize_count(interfaceDictionary):
    
    for pair in interfaceDictionary:
        
        for residue in interfaceDictionary[pair]:
            pdbCount = len(interfaceDictionary[pair][residue][2].split('$'))
            interfaceDictionary[pair][residue].append(interfaceDictionary[pair][residue][1] / pdbCount) #[3] is normalized by uniprot pair


# In[46]:


def average_histones(interfaceDictionary):
    
    avgDict = {}
    avgDict['type'] = {}
    
    for pair in interfaceDictionary:
        
        for residue in interfaceDictionary[pair]:                           
            targetFields = interfaceDictionary[pair][residue][0].split('@')[0].split('|')
            sourceFields = interfaceDictionary[pair][residue][0].split('@')[1].split('|') 

            if(targetFields[3] != 'other'): 
                histoneType = targetFields[3]
                normalizedCount = interfaceDictionary[pair][residue][3]
                
                if(histoneType in avgDict):
                    
                    if(residue in avgDict[histoneType]):
                        avgDict[histoneType][residue][0] += normalizedCount 
                        avgDict[histoneType][residue][1] += 1 #for normalization
                        
                    else:
                        avgDict[histoneType][residue] = [normalizedCount, 1]
                        
                else:
                    avgDict[histoneType] = {residue : [normalizedCount, 1]}
                    
            elif(sourceFields[3] != 'other'):
                histoneType = sourceFields[3]
                normalizedCount = interfaceDictionary[pair][residue][3]
                
                if(histoneType in avgDict):
                    
                    if(residue in avgDict[histoneType]):
                        avgDict[histoneType][residue][0] += normalizedCount 
                        avgDict[histoneType][residue][1] += 1 #for normalization
                        
                    else:
                        avgDict[histoneType][residue] = [normalizedCount, 1]
                        
                else:
                    avgDict[histoneType] = {residue : [normalizedCount, 1]}
    
    avgDict2 = {}
    avgDict2['type'] = {}
    
    for histoneType in avgDict:
        
        for residue in avgDict[histoneType]:
            
            if(histoneType in avgDict2):
                avgDict2[residue] = avgDict[histoneType][residue][0] / avgDict[histoneType][residue][1]

            else:
                avgDict2[histoneType] = {residue : avgDict[histoneType][residue][0] / avgDict[histoneType][residue][1]}

    return avgDict2


# In[48]:


def sum_contacts(interfaceDictionary):
    
    sumDict = {}
    histoneDict = {}
    
    for pair in interfaceDictionary:
        
        if(pair != 'uniprotPair'):
            
            targetFields = []
            sourceFields = []
            
            nucleosomeFlag = 0
            for residue in interfaceDictionary[pair]:
                pairs = interfaceDictionary[pair][residue][0].split('$')
                
                for instance in pairs:
                    targetFields = instance.split('@')[0].split('|')
                    sourceFields = instance.split('@')[1].split('|') 
                    
                    if(targetFields[7].split(':')[1] == '1' or sourceFields[7].split(':')[1] == '1'):
                        nucleosomeFlag = 1
                        break
                
                if(nucleosomeFlag):
                    break
            
            #randomResidue = list(interfaceDictionary[pair].keys())[0]
            #randomPair = interfaceDictionary[pair][randomResidue][0].split('$')[0]
            
            #targetFields = randomPair.split('@')[0].split('|')
            #sourceFields = randomPair.split('@')[1].split('|') 
            
            pdbList = []
            
            for residue in interfaceDictionary[pair]:
                pdbIDs = interfaceDictionary[pair][residue][2].split('$')
                
                for pdb in pdbIDs:
                    
                    if(pdb not in pdbList):
                        pdbList.append(pdb)
            
            if(targetFields[3] != 'other'):
                histoneType = targetFields[3]
                    
                newPair = histoneType + '@'
                
                if(sourceFields[3] != 'other'):
                    histoneType2 = sourceFields[3]
                                 
                    newPair += histoneType2
                
                else:   
                    
                    for field in sourceFields:
                        if(field == sourceFields[len(sourceFields) - 1]):
                            newPair += field
                        else:
                            newPair += field + '|'
                        
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount
                
                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]
                
                if(newPair in histoneDict):
                    histoneDict[newPair] += 1
                    
                else:
                    histoneDict = {newPair : 1}
                
            elif(sourceFields[3] != 'other'):
                histoneType = sourceFields[3]
                
                newPair = histoneType + '@'
                    
                for field in targetFields:
                    if(field == targetFields[len(targetFields) - 1]):
                        newPair += field
                    else:
                        newPair += field + '|'               
                    
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount

                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]
                  
                if(newPair in histoneDict):
                    histoneDict[newPair] += 1
                    
                else:
                    histoneDict = {newPair : 1}
                
            else:
                newPair = ''
                
                for field in targetFields:
                    if(field == targetFields[len(targetFields) - 1]):
                        newPair += field
                    else:
                        newPair += field + '|'
                    
                newPair += '@'
                
                for field in sourceFields:
                    if(field == sourceFields[len(sourceFields) - 1]):
                        newPair += field
                    else:
                        newPair += field + '|'  
                   
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount

                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]
                
                if(newPair in histoneDict):
                    histoneDict[newPair] += 1
                    
                else:
                    histoneDict = {newPair : 1}
                    
            
            
            
            
            for pair in histoneDict:
                sumDict[pair][0] = sumDict[pair][0] / histoneDict[pair]
            #####Have to account for the case when an interface between two non-histone chains is already in the dictionary, but is stored in a reverse order!!!!

    return sumDict


# In[49]:


def main():
    chainDictionary = {}
    chainDictionary['chain'] = {}
    get_chain_dictionaries(CHAIN_FILE, chainDictionary)
#     for pdb in chainDictionary:
#         for chain in chainDictionary[pdb]:
#             print(pdb + '\t' + chain + '\t' + str(chainDictionary[pdb][chain]))

    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')

    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    residue_count(interfaceFiles, chainDictionary, interfaceDictionary)
#     for pair in interfaceDictionary:
#         for a in interfaceDictionary[pair]:
#             print(pair + '\t' + a + '\t' + str(interfaceDictionary[pair][a]))

    normalize_count(interfaceDictionary)
#     for pair in interfaceDictionary:
#         for residue in interfaceDictionary[pair]:
#             uniprots = pair.split('@')
#             chains = str(interfaceDictionary[pair][residue][0]).split('@')
#             name1 = chains[0].split('|')[2]
#             name2 = chains[-1].split('|')[2]
#             count = str(interfaceDictionary[pair][residue][3])
#             print(uniprots[0] + '@' + name1 + '@' + uniprots[1] + '@' + name2 + '@' + residue + '@' + count)
    
#     avgDict = average_histones(interfaceDictionary)
#     for entry in avgDict:       
#         print(entry + '\t' + str(avgDict[entry]))
        
    sumDict = sum_contacts(interfaceDictionary)
    for pair in sumDict:
        chains = pair.split('@')
        target = chains[0]
        source = chains[1]
        
        pdbIDs = ''
        
        for pdb in sumDict[pair][1]:
            
            if(pdb == sumDict[pair][1][len(sumDict[pair][1]) - 1]):
                pdbIDs += pdb
                
            else:
                pdbIDs += pdb + '|'

            targetFields = target.split('|')
            sourceFields = source.split('|')      
            contacts = sumDict[pair][0]   
         ######
        if(len(targetFields) > 1 and targetFields[6].split(':')[1] == '1'):
            if(len(sourceFields) > 1):
                print(sourceFields[2] + ';' + sourceFields[0] + ';' + sourceFields[1] + ';' + sourceFields[4] + ';' + sourceFields[5] + ';' + sourceFields[6].split('nucleosome')[0] + ';' + targetFields[2] + ';' + targetFields[0] + ';' + targetFields[1] + ';' + targetFields[4] + ';' + targetFields[5] + ';' + targetFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'nucleosome' + ';' + str(contacts))
            else:
                print(source + ';' + ';' + ';' + ';' + ';' + ';' + targetFields[2] + ';' + targetFields[0] + ';' + targetFields[1] + ';' + targetFields[4] + ';' + targetFields[5] + ';' + targetFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'nucleosome' + ';' + str(contacts))   
        elif(len(sourceFields) > 1 and sourceFields[6].split(':')[1] == '1'):
            if(len(targetFields) > 1):
                print(sourceFields[2] + ';' + sourceFields[0] + ';' + sourceFields[1] + ';' + sourceFields[4] + ';' + sourceFields[5] + ';' + sourceFields[6].split('nucleosome')[0] + ';' + targetFields[2] + ';' + targetFields[0] + ';' + targetFields[1] + ';' + targetFields[4] + ';' + targetFields[5] + ';' + targetFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'nucleosome' + ';' + str(contacts))
            else:
                print(target + ';' + ';' + ';' + ';' + ';' + ';' + sourceFields[2] + ';' + sourceFields[0] + ';' + sourceFields[1] + ';' + sourceFields[4] + ';' + sourceFields[5] + ';' + sourceFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'nucleosome' + ';' + str(contacts))          
        elif(len(sourceFields) == 1 and len(targetFields) == 1):
            print(target + ';' + ';' + ';' + ';' + ';' + ';' + source + ';' + ';' + ';' + ';' + ';' + ';' + pdbIDs + ';' + 'NA' + ';' + str(contacts))
        else:
            if(len(targetFields) > 1 and len(sourceFields) > 1):
                print(sourceFields[2] + ';' + sourceFields[0] + ';' + sourceFields[1] + ';' + sourceFields[4] + ';' + sourceFields[5] + ';' + sourceFields[6].split('nucleosome')[0] + ';' + targetFields[2] + ';' + targetFields[0] + ';' + targetFields[1] + ';' + targetFields[4] + ';' + targetFields[5] + ';' + targetFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'histone' + ';' + str(contacts))
            elif(len(targetFields) > 1):
                print(source + ';' + ';' + ';' + ';' + ';' + ';' + targetFields[2] + ';' + targetFields[0] + ';' + targetFields[1] + ';' + targetFields[4] + ';' + targetFields[5] + ';' + targetFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'histone' + ';' + str(contacts))
            else:
                print(target + ';' + ';' + ';' + ';' + ';' + ';' + sourceFields[2] + ';' + sourceFields[0] + ';' + sourceFields[1] + ';' + sourceFields[4] + ';' + sourceFields[5] + ';' + sourceFields[6].split('nucleosome')[0] + ';' + pdbIDs + ';' + 'histone' + ';' + str(contacts))


# In[26]:


if __name__ == "__main__":
    main()

