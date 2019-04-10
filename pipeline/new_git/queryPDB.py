
# coding: utf-8

# In[61]:


import requests


# In[62]:


def main():
    url = "http://www.rcsb.org/pdb/rest/search"
    pdb_set = set()
    
    with open("combined/uniprots.txt", "r") as fh:
        for line in fh:
            uniprot = line.strip()
            query_text = """
<?xml version="1.0" encoding="UTF-8"?>

<orgPdbQuery>

<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>

<description>Simple query for a list of Uniprot Accession IDs:</description>

<accessionIdList>""" + uniprot + """</accessionIdList>

</orgPdbQuery>"""
            header = {"Content-Type": "application/x-www-form-urlencoded"}

            response = requests.post(url, data=query_text, headers=header)

            if response.status_code == 200:
                result = response.text.split('\n')
                
                #print("For %s found %d PDB entries matching query." % (uniprot, len(result) - 1))
                #print(response.text)
                for element in result:
                    if(element != ""):
                        pdb_set.add(element.split(':', 1)[0])
                    
            else:
                print("failed to retrieve results for %s" % uniprot)
                
    with open("combined/comboChains.txt", "w") as fh:                
        for pdb in pdb_set:
            query_text = "?pdbids=" + pdb + "&customReportColumns=compound,source,uniprotAcc,uniprotRecommendedName&service=wsdisplay&format=csv&ssa=null"
            urlCustom = "http://www.rcsb.org/pdb/rest/customReport"

            response = requests.post(urlCustom, data=query_text, headers=header)

            if(response.status_code == 200):
                    chains = response.text.split('<br />')
                    iterchains = iter(chains)
                    next(iterchains)

                    for chain in iterchains:
                        fh.write(chain + '\n')
                        
            else:
                print("failed to generate report for %s", pdb)


# In[ ]:


if __name__ == '__main__':
    main()

