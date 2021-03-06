
# coding: utf-8

# In[22]:


#!/usr/bin/env python 3
import re
import csv

#PATH = "./"
PATH = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Interfaces/"
CHAIN_FILE = "combined/comboChains.csv"
#CHAIN_FILE = "tempChains.csv"
PDB_LIST = "combined/comboID.txt"
#PDB_LIST = "tempID.txt"


# In[23]:


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


# In[24]:


#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone(name, typeCount):

    if(re.search(r'histone', name, re.I)):
       
        if(not re.search(r'chaperone|ase|ing|hsh49|rep|alpha|thioredoxin|bombinin|chromosomal|nucleoprotein|envelope|snrnp', name, re.I)):
            
            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain!!!)

            if(re.search(r'h2a', name, re.I)):
                typeCount[0] += 'H2A|'
                
            elif(re.search(r'h2b', name, re.I)):
                typeCount[0] += 'H2B|'
                    
            elif(re.search(r'h3', name, re.I)):
                typeCount[0] += 'H3|'
                    
            elif(re.search(r'h4', name, re.I)):
                typeCount[0] += 'H4|'
                    
            elif(re.search(r'h1', name, re.I)):
                typeCount[0] += 'H1|'  
                    
            elif(re.search(r'h5', name, re.I)):
                typeCount[0] += 'H5|'

            elif(re.search(r'arch', name, re.I)):
                typeCount[0] += 'archaeal histone|'
                
            else:
                typeCount[0] += 'some histone|'


# In[25]:


#PARAMETERS:
#name is a string with the name of chain
#typeCount is a 2-element list with the first one being an empty string and the second - 0 (iteger)

#RESULTS:
#Updates the typeCount array: the first element = general histone type(s); the second element = (should be) number of histones in chain

def is_histone2(name, typeCount):

    if(re.search(r'histone|h3|h4|h2a|lys\(ac\)', name, re.I)):
       
        if(not re.search(r'chaperone|ase|ing|fasciclin|rna|chain|region|fold|llama|anti|inhibitor|nanobody|nucleoprotein|regulatory|insulin|nr1h|upf0052|hiv|affitin|envelope|ddef|bh3|hsh49|1h4i|h2afy|h2a1', name, re.I)):
            
            typeCount[1] = 1 #adds the number of histones in chain  (should be changed to actual number of histones in chain!!!)

            if(re.search(r'h2a', name, re.I)):
                typeCount[0] += 'H2A|'
                
            elif(re.search(r'h2b', name, re.I)):
                typeCount[0] += 'H2B|'
                    
            elif(re.search(r'h4', name, re.I)):
                typeCount[0] += 'H4|'
                    
            elif(re.search(r'h1', name, re.I)):
                typeCount[0] += 'H1|'  
                    
            elif(re.search(r'h5', name, re.I)):
                typeCount[0] += 'H5|'

            elif(re.search(r'arch', name, re.I)):
                typeCount[0] += 'archaeal histone|'
                
            elif(re.search(r'h3|histone 3|', name, re.I)):
                typeCount[0] += 'H3|'  
                
            else:
                typeCount[0] += 'some histone|'


# In[26]:


#PARAMETERS: 
#pdbList is a text file with a header and one column PDB
#files is a list
#parameter is a string, either 'mapping' or 'interface' depending on desired results

#RESULTS:
#A list of absolute paths to either mapping files or interface files as it is stored on local NCBI machines


def get_files(pdbList, files, parameter):
    
    with open(pdbList, 'r') as pfh:
        
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


# In[27]:


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


# In[28]:


#PARAMETERS:
#cFile is tab-separated file with a header and 4 columns: pdb, chain, uniprot, name
#dictionary is nested with the innermost dict being dictionary['pdb'] = {}

#RESULTS: 
#Example: {1alq : {'G': 'A|p84233|Histone H3.2|H3|Xenopus laevis|nucleosome:1|bp:1|1|1|135'}}

def get_chain_dictionaries(cFile, dictionary): 
    
    #with open(cFile, 'r') as cfh:
        #cfh.readline()

    histoneCount = {} #is used to count number of histones in a structure!!!!!!!

    tempDict = {}
    tempDict['pdb'] = {}

    tempDict2 = {}
    tempDict2['pdb'] = {}
    
    mappingFiles = [] 
    get_files(PDB_LIST, mappingFiles, 'mapping')


    for file in mappingFiles:

        try: #adds a pdb entry to the dict only if mapping file exists

            with open(file, 'r') as mfh:
                mfh.readline() #skips header
                pdb = file.split('/')[-1].split('_', 1)[0]

                for mLine in mfh:
                    chainPair = mLine.split('\t')
                    alexChain = chainPair[0]
                    myChain = chainPair[1]
                    alignment = int(chainPair[7]) - int(chainPair[3])
                    mmCIFstart = int(chainPair[3])
                    mmCIFend = int(chainPair[5])
                    
                    if(pdb in tempDict):
                        if(alexChain in tempDict[pdb]):
                            tempDict2[pdb][alexChain].extend([alignment, mmCIFstart, mmCIFend])
                        else:
                            tempDict[pdb][alexChain] = myChain
                            tempDict2[pdb][alexChain] = [alignment, mmCIFstart, mmCIFend]
                        
                    else:
                        tempDict[pdb] = {alexChain : myChain}
                        tempDict2[pdb] = {alexChain : [alignment, mmCIFstart, mmCIFend]}

        except IOError:
            pass
            #print("Error: " + mappingFile + " does not appear to exist.")


    chains_file = csv.DictReader(open(cFile)) 

    for cLine in chains_file:

        pdb = cLine["structureId"].lower()
        if(pdb in tempDict): #continues only if a mapping file exists
            chain = cLine["chainId"]
            organism = cLine["source"]
            uniprot = cLine["uniprotAcc"].lower()
            name = cLine["uniprotRecommendedName"]

            histoneTypeAndCount = ['', 0]

            if(name):
                is_histone(name, histoneTypeAndCount) #checks whether the name looks like a histone!!!!!!!!!

            else:
                name = cLine["compound"]
                is_histone2(name, histoneTypeAndCount)

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
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType + organism + '|' #!!!!

                    else: #!!!!
                        dictionary[pdb][alexChain] = str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|' + organism  + '|' #!!!!

                else:

                    if(tempCount): #checks if chain is a histone!!!!
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + tempType + organism + '|'} #!!!!

                    else: #!!!!
                        dictionary[pdb] = {alexChain : str(tempDict[pdb][alexChain]) + '|' + uniprot + '|' + name + '|' + 'other|' + organism + '|'} #!!!!

            except ValueError:
                #print("Error: " + str(ValueError) + ", in " + pdb)               
                pass


    with open('combined/temp.tsv', 'w') as rfh:
        for structure in list(dictionary): 

            if(structure in histoneCount):
                uniqueHistoneNum = len(histoneCount[structure])
                partnerFlag = 0

                if(uniqueHistoneNum > 2): #checks if pdb has at least 3 different histones ~ is a nucleosome!!!!!!!

                    for chain in dictionary[structure]:
                        chainFields = dictionary[structure][chain].split('|')

                        chainType = chainFields[3]  ##

                        dictionary[structure][chain] += 'nucleosome:1|' #!!!!!!

                        if(partnerFlag == 0 and chainType == 'other'):
                            partnerFlag = 1

                    if(partnerFlag == 0):
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:0|'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '|'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 3)
                            rfh.write(str(dictionary[structure][chain]) + '\n')
                            
                            #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'no' + '\n')
                    else:
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1|'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '|'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 3)
                            rfh.write(str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'yes' + '\n')
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
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '|'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 3)
                            rfh.write(str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'histone' + '\t' + 'no' + '\n')
                    else:
                        for chain in dictionary[structure]:
                            dictionary[structure][chain] += 'bp:1|'
                            for element in tempDict2[structure][chain]:
                                dictionary[structure][chain] += str(element)
                                dictionary[structure][chain] += '|'
                            dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 3)
                            rfh.write(str(dictionary[structure][chain]) + '\n')
                            #rfh.write(structure + '\t' + 'histone' + '\t' + 'yes' + '\n')
            else:
                for chain in dictionary[structure]:
                    chainFields = dictionary[structure][chain].split('|')

                    chainType = chainFields[3]  ##

                    dictionary[structure][chain] += 'other:0|' #!!!!!!
                            
                    #rfh.write(structure + '\t' + 'nucleosome' + '\t' + 'no' + '\n')

                for chain in dictionary[structure]:
                    dictionary[structure][chain] += 'bp:1|'
                    for element in tempDict2[structure][chain]:
                        dictionary[structure][chain] += str(element)
                        dictionary[structure][chain] += '|'
                    dictionary[structure][chain] += str(len(tempDict2[structure][chain]) / 3)
                    rfh.write(str(dictionary[structure][chain]) + '\n')
                #del dictionary[structure]


# In[29]:


#PARAMETERS:
#interfaceFiles is a list of strings containing names of interface files
#chainDictionary is a dictionary produced by 'get_chain_dictionaries'
#interfaceDictonary is a nested dictionary with the innermost dict being interfaceDictionary['uniprotPair'] = {}, where the values are lists of the next form [chain pair data, residue number, pdb IDs] 
#Example: 'P84233@P62799': {'44': ['A|P84233|Histone H3.2|h3|1@B|P62799|Histone H4|h4|1$E|P84233|Histone H3.2|h3|1@F|P62799|Histone H4|h4|1', 17, '1zla$q5cl']
#note that data fields of a chain is separated by |, @ separates chain from binding partner chain, $ separates chain pairs and pdb structures
#A|P84233|Histone H3.2|H3|Xenopus laevis|nucleosome:1|bp:1


def residue_count(interfaceFiles, chainDictionary, interfaceDictionary):
    
    for file in interfaceFiles:
        pdb = file.split('/')[-1].split('_', 1)[0] 

        try:
            
            with open (file, 'r') as ifh:
                ifh.readline()

                with open('combined/results.tsv', 'a') as rfh:
                                     
                    for line in ifh:
                        lineFields = line.split('\t')

                        chain1 = lineFields[0].split('_', 1)[0] # the split part treats biological assembly chains as separate chains ???
                        chain2 = lineFields[4].split('_', 1)[0]

                        fields1 = chainDictionary[pdb][chain1].split('|')
                        fields2 = chainDictionary[pdb][chain2].split('|')

                        myChain1 = fields1[0]
                        myChain2 = fields2[0]
                        
                        type1 = fields1[3]
                        type2 = fields2[3]
                        
                        #nucleosome = int(fields1[5].split(':')[1])
                        bp = int(fields1[6].split(':')[1])
                        
                        align1 = int(fields1[7])
                        align2 = int(fields2[7])
                        
                        mmCIFstart1 = int(fields1[8])
                        mmCIFstart2 = int(fields2[8])
                        
                        mmCIFend1 = int(fields1[9])
                        mmCIFend2 = int(fields2[9])
                        
                        uniprotStart1 = mmCIFstart1 + align1
                        uniprotStart2 = mmCIFstart2 + align2
                        
                        uniprotEnd1 = mmCIFend1 + align1
                        uniprotEnd2 = mmCIFend2 + align2
                        
                        #if(not nucleosome and bp and type1 != 'other' and type2 == 'other'):
                        if(type1 != 'other' and type2 == 'other'):
                            
                            residue1 = lineFields[2]
                            residue2 = lineFields[6]
                            if((type1 == 'H1' and 1 <= uniprotStart1 <= 229 and 1 <= uniprotEnd1 <= 229) or (type1 == 'H2A' and 1 <= uniprotStart1 <= 385 and 1 <= uniprotEnd1 <= 385) or (type1 == 'H2B' and 1 <= uniprotStart1 <= 133 and 1 <= uniprotEnd1 <= 133) or (type1 == 'H3' and 1 <= uniprotStart1 <= 184 and 1 <= uniprotEnd1 <= 184) or (type1 == 'H4' and 1 <= uniprotStart1 <= 103 and 1 <= uniprotEnd1 <= 103)):
                                if(int(residue1) >= mmCIFstart1 and int(residue1) <= mmCIFend1):

                                    with open('combined/hitdata.txt', 'r') as hfh:
                                        hfh.readline()
                                        domainFlag = 0
                                        
                                        for line in hfh:
                                            fields = line.split('\t')


                                            pdb2 = fields[0].split(' - ')[1][0:4].lower() ###
                                            chain = fields[0].split(' - ')[1][4]  ###


                                            start = int(fields[3])
                                            end = int(fields[4])
                                            domtype = fields[1]
                                            name = fields[8]

                                            if(pdb2 == pdb and chain == myChain2):

                                                if(domtype == 'specific'):
                                                    if(int(residue2) >= int(start) and int(residue2) <= int(end)):
                                                        hfields = chainDictionary[pdb][chain1].split('|')
                                                        pfields = chainDictionary[pdb][chain2].split('|')
                                                        alignedRes = int(residue1) + align1
                                                        domainFlag = 1
                                                        rfh.write(hfields[1] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '_' + str(alignedRes) + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes) + '\t' + pdb + '\n')
                                                        #print(hfields[2] + '\t' + hfields[3] + '\t' + hfields[4] + '\t' + residue1 + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + pfields[4] + '\t' + name)
                                        hfh.seek(0)
                                        if(not domainFlag):
                                            hfields = chainDictionary[pdb][chain1].split('|')
                                            pfields = chainDictionary[pdb][chain2].split('|')
                                            alignedRes = int(residue1) + align1
                                            rfh.write(hfields[1] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '_' + str(alignedRes) + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes) + '\t' + pdb + '\n')
                                #else:
                                    #print(pdb, type1, residue1, mmCIFstart1, mmCIFend1)
                            #else:
                                #print(pdb + '\t' + fields1[0] + '\t' + fields1[2] + '\t' + str(uniprotStart1) + '\t' + str(uniprotEnd1) + '\n')

                        #if(not nucleosome and bp and type2 != 'other' and type1 == 'other'):
                        elif(type2 != 'other' and type1 == 'other'):

                            residue1 = lineFields[2]
                            residue2 = lineFields[6]
                            if((type2 == 'H1' and 1 <= uniprotStart2 <= 229 and 1 <= uniprotEnd2 <= 229) or (type2 == 'H2A' and 1 <= uniprotStart2 <= 385 and 1 <= uniprotEnd2 <= 385) or (type2 == 'H2B' and 1 <= uniprotStart2 <= 133 and 1 <= uniprotEnd2 <= 133) or (type2 == 'H3' and 1 <= uniprotStart2 <= 184 and 1 <= uniprotEnd2 <= 184) or (type2 == 'H4' and 1 <= uniprotStart2 <= 103 and 1 <= uniprotEnd2 <= 103)):
                                if(int(residue2) >= mmCIFstart2 and int(residue2) <= mmCIFend2):

                                    with open('combined/hitdata.txt', 'r') as hfh:
                                        hfh.readline()
                                        domainFlag = 0
                                        for line in hfh:
                                            fields = line.split('\t')


                                            pdb2 = fields[0].split(' - ')[1][0:4].lower() ###
                                            chain = fields[0].split(' - ')[1][4]  ###

                                            start = int(fields[3])
                                            end = int(fields[4])
                                            domtype = fields[1]
                                            name = fields[8]
                                            if(pdb2 == pdb and chain == myChain1):

                                                if(domtype == 'specific'):
                                                    if(int(residue1) >= int(start) and int(residue1) <= int(end)):
                                                        hfields = chainDictionary[pdb][chain2].split('|')
                                                        pfields = chainDictionary[pdb][chain1].split('|')
                                                        alignedRes = int(residue2) + align2
                                                        rfh.write(hfields[1] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '_' + str(alignedRes) + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + name + '\t' + pfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes) + '\t' + pdb + '\n')
                                                        domaingFlag = 1
                                                        #print(hfields[2] + '\t' + hfields[3] + '\t' + hfields[4] + '\t' + residue2 + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + pfields[4] + '\t' + name)
                                        hfh.seek(0)
                                        if(not domainFlag):
                                            hfields = chainDictionary[pdb][chain2].split('|')
                                            alignedRes = int(residue2) + align2
                                            pfields = chainDictionary[pdb][chain1].split('|')
                                            rfh.write(hfields[1] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '_' + str(alignedRes) + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + 'NA' + '\t' + pfields[4] + '\t' + hfields[3] + '\t' + str(alignedRes) + '\t' + pdb + '\n')
                                #else:
                                    #print(pdb, type2, residue2, mmCIFstart2, mmCIFend2)
                            #else:
                                #print(pdb + '\t' + fields2[0] + '\t' + fields2[2] + '\t' + str(uniprotStart2) + '\t' + str(uniprotEnd2) + '\n')
                                
                        elif(type2 == 'other' and type1 == 'other'):

                            hfields = chainDictionary[pdb][chain2].split('|')
                            pfields = chainDictionary[pdb][chain1].split('|')

                            rfh.write(hfields[1] +'\t' + hfields[2] + '\t' + hfields[4] + '\t' + hfields[3] + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + '' + '\t' + pfields[4] + '\t' + hfields[3] + '\t' + '' + '\t' + pdb + '\n')
                            #print(hfields[2] + '\t' + hfields[3] + '\t' + hfields[4] + '\t' + residue2 + '\t' + pfields[1] + '\t' + pfields[2] + '\t' + pfields[4] + '\t' + name)
                        

        except (IOError, KeyError) as e:
            #print("Error: " + file + " does not appear to exist.")
            print(e)
            pass


# In[30]:


def main():
    chainDictionary = {}
    chainDictionary['chain'] = {}
    get_chain_dictionaries(CHAIN_FILE, chainDictionary)
    #for pdb in chainDictionary:
     #   for chain in chainDictionary[pdb]:
      #      print(pdb + '\t' + chain + '\t' + str(chainDictionary[pdb][chain]))

    interfaceFiles = []
    get_files(PDB_LIST, interfaceFiles, 'interface')

    interfaceDictionary = {}
    interfaceDictionary['uniprotPair'] = {}
    residue_count(interfaceFiles, chainDictionary, interfaceDictionary)
#     for pair in interfaceDictionary:
#         for a in interfaceDictionary[pair]:
#             print(pair + '\t' + a + '\t' + str(interfaceDictionary[pair][a]))


# In[ ]:


if __name__ == "__main__":
    main()

