#!/usr/bin/env python3



def map_chains(mappingFile, chainMap):
	with open(mappingFile, 'r') as mf:
		for line in mf:
			chainMap[line.split()[0]] = line.split()[1]


def get_chain_lists(labeledChainsFile, pdbList, histoneList, partnerList):
		with open(labeledChainsFile, 'r') as lf:
			for pdbID in pdbList:
				for line in lf:
					if(line.startswith(pdbID)):
						histoneList[:] = line.split('\t')[1].split(',')
						partnerList[:] = line.split('\t')[2].split(',')
						lf.seek(0)
						break


def map_split_chains(labeledChainsFile, mappingFile, pdbList, histoneDict, partnerDict):
	histoneList = []
	partnerList = []
	get_chain_lists(labeledChainsFile, pdbList, histoneList, partnerList)
	print(histoneList)

	with open(mappingFile, 'r') as mf:
		for line in mf:
			if(line.split()[1] in histoneList):
				histoneDict[line.split()[0]] = line.split()[1]
			elif(line.split()[1] in partnerList):
				partnerDict[line.split()[0]] = line.split()[1]


#def get_frequency(file_, residueDict):
#	with open(file_, 'r') as residueFile:
#		for line in residueFile:




def main():
	mappingFile = "../data/Interfaces/6buz_chain_protein_mapping.tab"
	labeledChainsFile = "../data/labeled_chains.tsv"
	pdbIDs = ["4ZUX"]
	#For the following dicts:
	#key = original chain from Alex files; value = PDB chain from labeled_chains.tsv
	histones = {} 
	partners = {}

	map_split_chains(labeledChainsFile, mappingFile, pdbIDs, histones, partners)







if __name__ == "__main__":
	main()
