#!/usr/bin/env python
# coding: utf-8

# In[14]:


import scanpy as sc



import networkx as nx
import numpy as np

class DendrogramQuery:
    #object for querying node distances in dendogram
    #
    #dendrogramData should be output of sc.tl.dendrogram
    def __init__(self, dendrogramData):
        self.dendrogramData = dendrogramData
        
        self.cell2node = {}
        for id in dendrogramData["dendrogram_info"]["leaves"]:
            self.cell2node[dendrogramData["dendrogram_info"]["ivl"][id]] = id

        #create graph and compute shortest paths between all nodes
        G = nx.Graph()
        G.add_nodes_from(list(range(2+int(np.max(dendrogramData["linkage"][:,[0,1]].astype(int))))))
        joins = dendrogramData["linkage"][:,[0,1]].astype(int)
        numInitalNodes = np.max(dendrogramData["dendrogram_info"]["leaves"])
        for i in range(joins.shape[0]):
            G.add_edge(joins[i][0],numInitalNodes+1+i)
            G.add_edge(joins[i][1],numInitalNodes+1+i)
            
        self.shortestPaths = dict(nx.all_pairs_shortest_path_length(G))
        self.G = G
    ### search should be done on values of group by
    def dist(self,cell1,cell2):
        return self.shortestPaths[self.cell2node[cell1]][self.cell2node[cell2]]




import pandas as pd
import numpy as np

class CellOntologyQuery:
    #object for querying node distances in dendogram
    #
    #dendrogramData should be output of sc.tl.dendrogram
    def parse_group(lines):
        ## added a small change to map all names to lower case
        term_id = lines[[o.startswith("id:") for o in lines]][0].strip().split()[1].lower()
        try:
            #print("-------------------")
            is_a = [s.strip().split()[1].lower() for s in lines[[o.startswith("is_a:") for o in lines]]]
            
            #is_a = lines[[o.startswith("is_a:") for o in lines]][0].strip().split()[1].lower()
        except:
            is_a = None
        name = lines[[o.startswith("name:") for o in lines]][0].strip().split("name:")[1].strip().lower()
        return name, is_a, term_id

    def __init__(self):
        #### read and parse cell ontology data into dataframe, parsing was borrowed from yanay
        obo_loc = "/dfs/project/cross-species/yanay/data/tabula/cl.obo.txt"
        with open(obo_loc, "r", encoding='utf-8') as f:
            obo = f.readlines()
        obo_term_idxs = np.where([o.startswith('[Term]') for o in obo])[0]

        
        all_rows = []
        for i in range(1, len(obo_term_idxs) - 1):
            ls = np.array(obo[obo_term_idxs[i]:obo_term_idxs[i+1]])
            r = CellOntologyQuery.parse_group(ls)
            if r[1] is None:
                all_rows.append(r)
            else:
                for parent in r[1]:
                    all_rows.append((r[0],parent,r[2]))
        
        obo_tbl = pd.DataFrame(all_rows, columns=["name", "is_a", "id"]).set_index("id")
        obo_tbl.index = obo_tbl.index.astype(str)
        #print(obo_tbl.groupby(["id"]).count().max())
        
        #print(obo_tbl.query("id == 'cl:0000235'"))
        
        G = nx.Graph()
        G.add_nodes_from(obo_tbl.index)
        
        # add edges to graph, also compute mapping between nodes and names of cells
        self.cell2node = {}
        for el in obo_tbl.iterrows():
            if el[1]["is_a"] is not None:
                G.add_edge(el[0],el[1]["is_a"])
                self.cell2node[el[1]["name"]] = el[0]
            else:
                 self.cell2node[el[1]["name"]] = el[0]
           
                
        
        self.shortestPaths = dict(nx.all_pairs_shortest_path_length(G))
        self.G = G
        
        self.notFound = set()
    ### search should be done on values of group by
    def dist(self,cell1,cell2):
        try:
            return self.shortestPaths[self.cell2node[cell1]][self.cell2node[cell2]]
        except:
            try:
                self.cell2node[cell1]
            except:
                self.notFound.add(cell1)#print(f"!!! did not fint key '{cell1}'")
            try:
                self.cell2node[cell2]
            except:
                self.notFound.add(cell2)#print(f"!!! did not fint key '{cell2}'")
                
            return np.inf


