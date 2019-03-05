#!/usr/bin/env python
# coding: utf-8

# In[45]:


#!/usr/bin/env python 3
import re
import csv
import collections

PATH = "../data/Interfaces/"
#PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "text.tsv"
PDB_LIST = "pdbList.txt"


# In[33]:


def is_histone(name):
   
    if(re.search(r'histone', name, re.I)):
       
        if(not re.search(r'chaperone|ase|ing|hsh49|rep|alpha|thioredoxin|bombinin|chromosomal|nucleoprotein|envelope|snrnp', name, re.I)): 
            
            if(re.search(r'h2a', name, re.I)):
                return'H2A'
                
            elif(re.search(r'h2b', name, re.I)):
                return 'H2B'
                    
            elif(re.search(r'h3', name, re.I)):
                return 'H3'
                    
            elif(re.search(r'h4', name, re.I)):
                return 'H4'
                    
            elif(re.search(r'h1', name, re.I)):
                return 'H1'  
                    
            elif(re.search(r'h5', name, re.I)):
                return 'H5'
             
            elif(re.search(r'arch', name, re.I)):
                return 'archaeal histone'
            
            else:
                return 'some histone'
            
    return 'binding partner'


# In[ ]:


#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS: 
#The format of the end-product dictionary is: {pdb : {AlexChain: myChain|UNIPROT|name|type|process|function|organism|nucleosome(bool)|bindingPartner(bool)}}
#Example: {1alq : {'G': 'E|P02302|Histone H3.3C|H3|1212|4141|human|1|0'}}


def get_chain_dictionaries(cFile, dictionary): 
    
    chains_file = csv.DictReader(open(cFile))
        
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


# In[42]:


def main():

    chains_file = csv.DictReader(open("chains.csv"))
    contain_histones = defaultdict(dict)
    
    for row in chains_file:
        
        if(rowtempCount): #if the chain is a [part of a] histone

            if(pdb in histoneCount):

                if(tempType not in histoneCount[pdb]):
                    histoneCount[pdb].append(tempType) #!!!!!!

            else:
                histoneCount[pdb] = [tempType]       




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


# In[43]:


if __name__ == "__main__":
    main()


# In[ ]:


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


# In[ ]:




