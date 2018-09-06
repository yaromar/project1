
# coding: utf-8

# In[21]:


def main():
    
    
    with open ('processedPDBlist.tsv', 'r') as pfh:
        pfh.readline() #skips header
        
        depDict = {}
        depDict['type'] = {}
        
        for line in pfh:
            lineFields = line.split('\t') #gets only the first 8 columns !!!
            year = lineFields[2].split('-', 1)[0]
            pdbType = lineFields[3].strip() + '#' + lineFields[4].strip()
            
            if(pdbType in depDict):
                
                if(year in depDict[pdbType]):
                    depDict[pdbType][year] += 1
                
                else:
                    depDict[pdbType][year] = 1
                    
            else:
                depDict[pdbType] = {year : 1}
                
        for pdbType in depDict:
            print('###################' + pdbType + '###################')
            
            for year in depDict[pdbType]:
                print(year + '\t' + str(depDict[pdbType][year]))


# In[22]:


if __name__ == "__main__":
    main()

