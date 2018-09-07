
# coding: utf-8

# In[17]:


def main():
    
    with open ('processedPDBlist.tsv', 'r') as pfh:
        pfh.readline() #skips header
        
        depDict = {}
        depDict['type'] = {}
        
        for line in pfh:
            lineFields = line.split('\t') #gets only the first 8 columns !!!
            year = lineFields[3].split('-', 1)[0]
            pdbType = lineFields[1].strip() + '#' + lineFields[2].strip()

            if(pdbType in depDict):
                
                if(year in depDict[pdbType]):
                    depDict[pdbType][year] += 1
                
                else:
                    depDict[pdbType][year] = 1
                    
            else:
                depDict[pdbType] = {year : 1}
                
        for pdbType in depDict:
            print('###################' + pdbType + '###################')
            cumSum = 0
            
            for year in depDict[pdbType]:
                cumSum += depDict[pdbType][year]
                print(year + '\t' + str(depDict[pdbType][year]) + '\t' + str(cumSum))
                


# In[18]:


if __name__ == "__main__":
    main()

