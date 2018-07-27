#!/usr/bin/env python 3
import os

PATH = "../data/Interfaces/"




def get_chain_lists(labeledChainsFile, mappedChains):
        with open(labeledChainsFile, 'r') as lfh:
            lfh.readline() #skips header

            for line in lfh:
                pdbID = line.split('\t', 1)[0]
                histoneChains = line.split('\t')[1].split(',')
                partnerChains = line.split('\t')[2].split(',')

                histDict = {}
                partnDict = {}
                for chain in histoneChains:
                    histDict[chain] = ''

                for chain in partnerChains:
                    partnDict[chain] = ''

                mappedChains[pdbID] = {'histone' : {}}
                mappedChains[pdbID] = {'partner' : {}}

                mappedChains[pdbID]['histone'] = histDict
                mappedChains[pdbID]['partner'] = partnDict


def map_chains(labeledChainsFile, mappingFiles, mappedChains):
    get_chain_lists(labeledChainsFile, mappedChains)
    
    for file in mappingFiles:
        pdbID = file.split('_', 1)[0]
        with open(PATH+file, 'r') as mfh:
            mfh.readline() #skips header  
            
            for line in mfh:
                lineFields = line.split('\t', 2)
                chainOriginal = lineFields[0] #Alexander's files
                chainNew = lineFields[1] #labeled_chains file
                
                if(chainNew in mappedChains[pdbID]['histone']):
                    mappedChains[pdbID]['histone'][chainNew] = chainOriginal
                elif(chainNew in mappedChains[pdbID]['partner']):
                    mappedChains[pdbID]['partner'][chainNew] = chainOriginal



class pdbFreq:
    def __init__(self, interfaceFiles, mappedChains):
        self.freq = {}
        self.freq['pdb'] = {}
        self.freq['pdb']['chain'] = {}
        self.freq['pdb']['chain']['residue'] = {}

        for file in interfaceFiles:
            pdbID = file.split('_', 1)[0]

            with open (PATH+file, 'r') as ifh:
                for line in ifh:
                    lineFields = line.split('\t', 7) #gets only the first 8 columns !!!
                    chain1 = lineFields[0]
                    chain2 = lineFields[4]

                    if((chain1 in mappedChains[pdbID]['histone'].values()) and (chain2 in mappedChains[pdbID]['partner'].values())):
                        res = lineFields[2]
                        self.addResidue(pdbID, chain1, res)
                    elif((chain1 in mappedChains[pdbID]['partner'].values()) and (chain2 in mappedChains[pdbID]['histone'].values())):
                        res = lineFields[6]
                        self.addResidue(pdbID, chain2, res)

    def addResidue(self, pdb, ch, aa):
        if(pdb in self.freq):
            if(ch in self.freq[pdb]):
                if(aa in self.freq[pdb][ch]):
                    self.freq[pdb][ch][aa] += 1
                else:
                    self.freq[pdb][ch][aa] = 1
            else:
                self.freq[pdb][ch] = {aa : 1}
        else:
            self.freq[pdb] = {ch : {aa : 1}}

    def printContent(self):
        print(self.freq)





def main():
    folder = os.listdir(PATH)
    
    mappingFiles = []
    interfaceFiles = []
    for file in folder:
        if("mapping" in file):
            mappingFiles.append(file)
        elif("contacts"in file):
            interfaceFiles.append(file)
    
    labeledChainsFiles = "../data/labeled_chains.tsv"

    mappedChains = {} 
    mappedChains['PDB'] = {}
    mappedChains['PDB']['type'] = {}
    mappedChains['PDB']['type']['chain'] = {}
    
    map_chains(labeledChainsFiles, mappingFiles, mappedChains)
    
    result = pdbFreq(interfaceFiles, mappedChains)
    result.printContent()




if __name__ == "__main__":
main()
