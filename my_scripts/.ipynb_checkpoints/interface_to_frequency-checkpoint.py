{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 231,
   "metadata": {},
   "outputs": [],
   "source": [
    "#!/usr/bin/env python 3"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 232,
   "metadata": {},
   "outputs": [],
   "source": [
    "#def map_chains(mappingFile, chainMap):\n",
    "#    with open(mappingFile, 'r') as mfh:\n",
    "#        for line in mfh:\n",
    "#            chainMap[line.split()[0]] = line.split()[1]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 297,
   "metadata": {},
   "outputs": [],
   "source": [
    "def get_chain_lists(labeledChainsFile, mappedChains):\n",
    "        with open(labeledChainsFile, 'r') as lfh:\n",
    "            lfh.readline() #skips header\n",
    "            \n",
    "            for line in lfh:\n",
    "                pdbID = line.split('\\t', 1)[0]\n",
    "                histoneChains = line.split('\\t')[1].split(',')\n",
    "                partnerChains = line.split('\\t')[2].split(',')\n",
    "                \n",
    "                histDict = {}\n",
    "                partnDict = {}\n",
    "                for chain in histoneChains:\n",
    "                    histDict[chain] = ''\n",
    "\n",
    "                for chain in partnerChains:\n",
    "                    partnDict[chain] = ''\n",
    "\n",
    "                mappedChains[pdbID] = {'histone' : {}}\n",
    "                mappedChains[pdbID] = {'partner' : {}}\n",
    "                \n",
    "                mappedChains[pdbID]['histone'] = histDict\n",
    "                mappedChains[pdbID]['partner'] = partnDict"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 298,
   "metadata": {},
   "outputs": [],
   "source": [
    "def map_chains(labeledChainsFile, mappingFiles, mappedChains):\n",
    "    get_chain_lists(labeledChainsFile, mappedChains)\n",
    "    \n",
    "    for file in mappingFiles:\n",
    "        pdbID = file.split('/', 3)[3]\n",
    "        pdbID = pdbID.split('_', 1)[0]\n",
    "        \n",
    "        with open(file, 'r') as mfh:\n",
    "            mfh.readline() #skips header  \n",
    "            \n",
    "            for line in mfh:\n",
    "                lineFields = line.split('\\t', 2)\n",
    "                chainOriginal = lineFields[0] #Alexander's files\n",
    "                chainNew = lineFields[1] #labeled_chains file\n",
    "                \n",
    "                if(chainNew in mappedChains[pdbID]['histone']):\n",
    "                    mappedChains[pdbID]['histone'][chainNew] = line.split()[0]\n",
    "                elif(chainNew in mappedChains[pdbID]['partner']):\n",
    "                    mappedChains[pdbID]['partner'][chainNew] = line.split()[0]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 299,
   "metadata": {},
   "outputs": [],
   "source": [
    "class pdbFreq:\n",
    "    def __init__(self, interfaceFiles, mappedChains):\n",
    "        self.freq = {}\n",
    "        self.freq['chain'] = {}\n",
    "        self.freq['chain']['residue'] = {}\n",
    "            \n",
    "        for file in interfaceFiles:\n",
    "            pdbID = file.split('/', 3)[3]\n",
    "            pdbID = pdbID.split('_', 1)[0]\n",
    "            \n",
    "            with open (file, 'r') as ifh:\n",
    "                for line in ifh:\n",
    "                    lineFields = line.split('\\t', 7) #gets only the first 8 columns !!!\n",
    "                    chain1 = lineFields[0]\n",
    "                    chain2 = lineFields[4]\n",
    "                \n",
    "                    if((chain1 in mappedChains[pdbID]['histone'].values()) and (chain2 in mappedChains[pdbID]['partner'].values())):\n",
    "                        res = lineFields[2]\n",
    "                        self.addResidue(chain1, res)\n",
    "                    elif((chain1 in mappedChains[pdbID]['partner'].values()) and (chain2 in mappedChains[pdbID]['histone'].values())):\n",
    "                        res = lineFields[6]\n",
    "                        self.addResidue(chain2, res)\n",
    "                    \n",
    "    def addResidue(self, ch, aa):\n",
    "        if(ch in self.freq):\n",
    "            if(aa in self.freq[ch]):\n",
    "                self.freq[ch][aa] += 1\n",
    "            else:\n",
    "                self.freq[ch][aa] = 1\n",
    "        else:\n",
    "            self.freq[ch] = {aa: 1}\n",
    "            \n",
    "    def printContent(self):\n",
    "        print(self.freq)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 300,
   "metadata": {
    "scrolled": true
   },
   "outputs": [],
   "source": [
    "def main():\n",
    "    labeledChainsFiles = \"../data/labeled_chains.tsv\"\n",
    "    mappingFiles = [\"../data/Interfaces/4zux_chain_protein_mapping.tab\"]\n",
    "    interfaceFiles = [\"../data/Interfaces/4zux_atomic_contacts_5.0A.tab\"]\n",
    "\n",
    "    mappedChains = {} \n",
    "    mappedChains['PDB'] = {}\n",
    "    mappedChains['PDB']['type'] = {}\n",
    "    mappedChains['PDB']['type']['chain'] = {}\n",
    "    \n",
    "    map_chains(labeledChainsFiles, mappingFiles, mappedChains)\n",
    "    \n",
    "    result = pdbFreq(interfaceFiles, mappedChains)\n",
    "    result.printContent()\n",
    "    #test = populate_pdb_freqs(interfaceFiles, histones, partners)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 301,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'chain': {'residue': {}}, 'A': {'81': 2}, 'B': {'80': 11}, 'G': {'25': 19, '23': 4, '57': 7, '72': 27, '65': 20, '69': 17, '68': 5, '61': 8, '58': 26, '62': 30}, 'H': {'45': 56, '47': 22, '46': 10, '48': 15, '49': 1, '114': 21, '111': 26, '107': 22, '118': 18, '117': 36, '110': 8, '113': 14, '120': 3, '121': 18}, 'C': {'58': 28, '20': 2, '23': 9, '25': 24, '57': 2, '61': 9, '62': 29, '65': 23, '68': 2, '69': 20, '72': 19}, 'D': {'107': 27, '108': 1, '110': 19, '111': 32, '113': 17, '114': 26, '116': 12, '117': 44, '118': 21, '120': 4, '121': 7, '44': 1, '45': 62, '46': 7, '47': 21, '48': 14}}\n"
     ]
    }
   ],
   "source": [
    "if __name__ == \"__main__\":\n",
    "    main()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
