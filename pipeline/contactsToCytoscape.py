
# coding: utf-8

# In[42]:


#NOTES:
#To get a list of PDB that contains histones from 'text.csv':
#cut -f1 text.tsv | uniq | awk '{print tolower($0)}' | sort

#text.tsv should be sorted by pdb and then uniprot name
#it MUST have 'NA' in uniprot and name blanks


# In[43]:


#!/usr/bin/env python 3
import re

#PATH = "../data/Interfaces/"
PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "text.tsv"
PDB_LIST = "pdbList.txt"


# In[44]:


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


# In[45]:


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
                typeCount[0] += 'some histone#'


# In[46]:


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


# In[47]:


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


# In[97]:


#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS:
#The format of the end-product dictionary is: {pdb : {AlexChain: myChain#UNIPROT#name#type#nucleosome(bool)#bindingPartner(bool)}}
#Example: {1alq : {'G': 'E#P02302#Histone H3.3C#H3#1#0'}}


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
                
                
                #######################
                if(tempCount):
                    
                    if(pdb in histoneCount):
                        
                        if(tempType not in histoneCount[pdb]):
                            histoneCount[pdb].append(tempType) #!!!!!!

                    else:
                        histoneCount[pdb] = [tempType]

                        
                        
                        
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
            partnerFlag = 0 ###
            uniqueHistoneNum = len(histoneCount[structure])
            
            if(uniqueHistoneNum > 3): #checks if pdb has at least a half of nucleosome!!!!!!!
                
                for chain in dictionary[structure]:
                    chainFields = dictionary[structure][chain].split('#')

                    chainType = chainFields[-2]
                    chainName = chainFields[1] + '\t' + chainFields[2]

                    dictionary[structure][chain] += '1#' #!!!!!!
                    
                    if(partnerFlag == 0 and chainType == 'other'):
                        #print('\t\t\t\t\t' + structure + '\t' + 'nucleosome' + '\t' + chainName)
                        partnerFlag = 1
                    
                if(partnerFlag == 0):
                    dictionary[structure][chain] += '0'
                    print(structure + '\t' + 'nucleosome' + '\t' + 'no bp')
                    
                else:
                    dictionary[structure][chain] += '1'
                    print(structure + '\t' + 'nucleosome' + '\t' + 'yes bp')
            
            else: #!!!!!!

                for chain in dictionary[structure]: #!!!!!
                    chainFields = dictionary[structure][chain].split('#')
                    
                    chainType = chainFields[-2]
                    chainName = chainFields[1] + '\t' + chainFields[2]
                    
                    dictionary[structure][chain] += '0#' #!!!!!! 

                    if(partnerFlag == 0 and chainType == 'other'):
                        #print('\t\t\t\t\t' + structure + '\t' + 'histone' + '\t' + chainName)
                        partnerFlag = 1

                if(partnerFlag == 0):
                    dictionary[structure][chain] += '0'
                    print(structure + '\t' + 'histone' + '\t' + 'no bp')

                else:
                    dictionary[structure][chain] += '1'
                    print(structure + '\t' + 'histone' + '\t' + 'yes bp')


# In[98]:


#PARAMETERS:
#interfaceFiles is a list of strings containing names of interface files
#chainDictionary is a dictionary produced by 'get_chain_dictionaries'
#interfaceDictonary is a nested dictionary with the innermost dict being interfaceDictionary['uniprotPair'] = {}, where the values are lists of the next form [chain pair data, residue number, pdb IDs] 
#Example: 'P84233@P62799': {'44': ['A#P84233#Histone H3.2#h3#1@B#P62799#Histone H4#h4#1$E#P84233#Histone H3.2#h3#1@F#P62799#Histone H4#h4#1', 17, '1zla$q5cl']
#note that data fields of a chain is separated by #, @ separates chains binding partner chains, $ separates chain pairs and pdb structures


def residue_count(interfaceFiles, chainDictionary, interfaceDictionary):
    for file in interfaceFiles:
        pdb = file.split('/')[-1].split('_', 1)[0] 
        
        try:
            with open (file, 'r') as ifh:
                ifh.readline() #skips header

                for line in ifh:
                    lineFields = line.split('\t', 6) #gets only the first 8 columns !!!

                    chain1 = lineFields[0].split('_', 1)[0] # the split part treats biological assembly chains as separate chains ???
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
                                
                                if(pdb not in interfaceDictionary[uniprotPair1][residue1][2]):
                                    interfaceDictionary[uniprotPair1][residue1][2] += '$' + pdb
                                   
                        else:
                            interfaceDictionary[uniprotPair1][residue1] = [chainPair1, 1, pdb]    #format of the innermost dict        

                    elif(uniprotPair2 in interfaceDictionary):

                        if(residue2 in interfaceDictionary[uniprotPair2]):
                            interfaceDictionary[uniprotPair2][residue2][1] += 1

                            if(chainPair2 not in interfaceDictionary[uniprotPair2][residue2][0]):
                                interfaceDictionary[uniprotPair2][residue2][0] += '$' + chainPair2                        

                                if(pdb not in interfaceDictionary[uniprotPair2][residue2][2]):
                                    interfaceDictionary[uniprotPair2][residue2][2] += '$' + pdb
                                    
                        else:
                            interfaceDictionary[uniprotPair2][residue2] = [chainPair2, 1, pdb]

                    else:
                         interfaceDictionary[uniprotPair1] = {residue1 : [chainPair1, 1, pdb]}  
                            
        except IOError:
            #print("Error: " + interfaceFile + " does not appear to exist.")
            pass


# In[99]:


def normalize_count(interfaceDictionary):
    
    for pair in interfaceDictionary:
        
        for residue in interfaceDictionary[pair]:
            pdbCount = len(interfaceDictionary[pair][residue][2].split('$'))
            interfaceDictionary[pair][residue].append(interfaceDictionary[pair][residue][1] / pdbCount) #[3] is normalized by uniprot pair


# In[100]:


def average_histones(interfaceDictionary):
    
    avgDict = {}
    avgDict['type'] = {}

    for pair in interfaceDictionary:
        
        for residue in interfaceDictionary[pair]:                           
            targetFields = interfaceDictionary[pair][residue][0].split('@')[0].split('#')
            sourceFields = interfaceDictionary[pair][residue][0].split('@')[-1].split('#')
            
            if(targetFields[-2] != 'other' and targetFields[-3] != 'other'): #MAKE ENTRIES HAVE THE SAME NUMBER OF ELEMENTS!!!
                histoneType = targetFields[3]
                normalizedCount = interfaceDictionary[pair][residue][3]
                
                if(histoneType in avgDict):
                    
                    if(residue in avgDict[histoneType]):
                        avgDict[histoneType][residue] += normalizedCount
                    
                    else:
                        avgDict[histoneType][residue] = normalizedCount
                        
                else:
                    avgDict[histoneType] = {residue : normalizedCount}
                    
            elif(sourceFields[-2] != 'other' and sourceFields[-3] != 'other'):#MAKE ENTRIES HAVE THE SAME NUMBER OF ELEMENTS!!!
                histoneType = sourceFields[3]
                normalizedCount = interfaceDictionary[pair][residue][3]
                
                if(histoneType in avgDict):
                    
                    if(residue in avgDict[histoneType]):
                        avgDict[histoneType][residue] += normalizedCount
                    
                    else:
                        avgDict[histoneType][residue] = normalizedCount
                        
                else:
                    avgDict[histoneType] = {residue : normalizedCount}
    
    return avgDict


# In[1]:


def sum_contacts(interfaceDictionary):
    
    sumDict = {}

    for pair in interfaceDictionary:
        if(pair != 'uniprotPair'):

            randomResidue = list(interfaceDictionary[pair].keys())[0]

            targetFields = interfaceDictionary[pair][randomResidue][0].split('@')[0].split('#')
            sourceFields = interfaceDictionary[pair][randomResidue][0].split('@')[-1].split('#')
            
            pdbList = []
            for residue in interfaceDictionary[pair]:
                pdbIDs = interfaceDictionary[pair][residue][2].split('$')
                
                for pdb in pdbIDs:
                    
                    if(pdb not in pdbList):
                        pdbList.append(pdb)
            
            #print(pair, '\t', pdbList)
            
            interfaceFlag = 0
            
            if(targetFields[-2] != 'other' and targetFields[-3] != 'other'):
                histoneType = targetFields[3]
                    
                newPair = histoneType + '@'
                
                if(sourceFields[-2] != 'other' and sourceFields[-3] != 'other'):
                    histoneType2 = sourceFields[3]
                                 
                    newPair += histoneType2
                    
                else:   
                    for field in sourceFields:
                        newPair += field + '#'
                        interfaceFlag = 1
                        
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount
                
                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]

            elif(sourceFields[-2] != 'other' and sourceFields[-3] != 'other'):
                histoneType = sourceFields[3]
                
                newPair = histoneType + '@'
                    
                for field in targetFields:
                    newPair += field + '#'                 
                    
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount

                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]
                    
            else:
                newPair = ''
                
                for field in targetFields:
                    newPair += field + '#'
                    
                newPair += '@'
                
                for field in sourceFields:
                    newPair += field + '#'    
                   
                totalCount = 0

                for residue in interfaceDictionary[pair]:
                    normalizedCount = interfaceDictionary[pair][residue][3]
                    totalCount += normalizedCount

                if(newPair in sumDict):
                    sumDict[newPair][0] += totalCount
                
                else:
                    sumDict[newPair] = [totalCount, pdbList]
            
            #####Have to account for the case when an interface between two non-histone chains is already in the dictionary, but is stored in a reverse order!!!!

    return sumDict


# In[102]:


def main():
    chainDictionary = {}
    chainDictionary['chain'] = {}
    
    get_chain_dictionaries(CHAIN_FILE, chainDictionary)
    
    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')
    
    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    
    residue_count(interfaceFiles, chainDictionary, interfaceDictionary)
    
    normalize_count(interfaceDictionary)
    #for pair in interfaceDictionary:
    #    print(pair + ': ' + str(interfaceDictionary[pair]))
    #    print('\n')
    
    avgDict = average_histones(interfaceDictionary)
    #print(avgDict)    
    
    sumDict = sum_contacts(interfaceDictionary)
    print('target' + '\t' + 'source' + '\t' + 'contacts')
    pdbList = '1zla, 1aoi,1eqz,1f66,1hio,1hq3,1id3,1kx3,1kx4,1kx5,1m18,1m19,1m1a,1p34,1p3a,1p3b,1p3f,1p3g,1p3i,1p3k,1p3l,1p3m,1p3o,1p3p,1s32,1tzy,1u35,1zbb,2aro,2cv5,2f8n,2fj7,2hio,2nqb,2nzd,2pyo,3a6n,3afa,3an2,3av1,3av2,3ayw,3aze,3azf,3azg,3azh,3azi,3azj,3azk,3azl,3azm,3azn,3b6f,3b6g,3c1b,3c1c,3kuy,3kwq,3kxb,3lel,3lja,3lz0,3lz1,3mgp,3mgq,3mgr,3mgs,3mnn,3mvd,3o62,3reh,3rei,3rej,3rek,3rel,3tu4,3ut9,3uta,3utb,3w96,3w97,3w98,3w99,3wa9,3waa,3wkj,3wtp,3x1s,3x1t,3x1u,3x1v,4j8u,4j8v,4j8w,4j8x,4jjn,4kgc,4kud,4ld9,4qlc,4r8p,4wu8,4wu9,4x23,4xuj,4xzq,4ym5,4ym6,4ys3,4z5t,4z66,4zux,5av5,5av6,5av8,5av9,5avb,5avc,5ay8,5b0y,5b0z,5b1l,5b1m,5b24,5b2i,5b2j,5b31,5b32,5b33,5b40,5cp6,5cpi,5cpj,5cpk,5dnm,5dnn,5e5a,5f99,5gse,5gsu,5gt0,5gt3,5gtc,5gxq,5hq2,5jrg,5kgf,5mlu,5nl0,5o9g,5omx,5ong,5onw,5oxv,5oy7,5x0x,5x0y,5x7x,5xf3,5xf4,5xf5,5xf6,5xm0,5xm1,6buz,6c0w,6esf,6esg,6esh,6esi,6etx,6fml,6fq5,6fq6,6fq8'.split(',')

    for pair in sumDict:
        chains = pair.split('@')
        target = chains[0]
        source = chains[1]
        
        contacts = sumDict[pair][0]
        pdbIDs = ''
        
        for pdb in sumDict[pair][1]:
            pdbIDs += pdb + '#'

        targetFields = target.split('#')
        sourceFields = source.split('#')
        
        smallList = pdbIDs.split('#')
        
        intersectionFlag = 0
        
        for sPDB in smallList:
            
            if(sPDB in pdbList):
                intersectionFlag = 1
                break
        
        if(len(targetFields) > 1 and intersectionFlag):
            print(source + '\t' + target + '\t' + pdbIDs + '\t' + 'nucleosome')
            
        elif(len(sourceFields) > 1 and intersectionFlag):
            print(target + '\t' + source + '\t' + pdbIDs + '\t' + 'nucleosome')
            
        elif(intersectionFlag):
            print(target + '\t' + source + '\t' + pdbIDs + '\t' + 'nucleosome')
            
        else:
            print(target + '\t' + source + '\t' + pdbIDs + '\t' + 'histone')
            
        #print(target + '\t' + source + '\t' + str(contacts))
        


# In[103]:


if __name__ == "__main__":
    main()

