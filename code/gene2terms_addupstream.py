import networkx as nx
from collections import defaultdict
import pymysql 



#############################################################################################################
#
#    SQL STATEMENT EXECUTED IN PYTHON
#    RETURNS THE PPI AS A NETWORKX GRAPH 
#    OPTIONS:
#            - G = get_ppi(0)  CONTAINS ALL CONNECTED COMPONENTS
#            - G = get_ppi(1)  CONTAINS ONLY LCC (Largest Connected Component)
#
#    IMPLEMENTED FILTERS ARE:
#            - NO SELFLOOPS
#            - ONLY CONNECTIONS THAT HAVE AN AUTHOR AND PubMed ID ARE INCLUDED
#            - ONLY PROTEIN CODING GENES ARE INCLUDED 
#      (CONTAINS Locus Types: 'gene with protein product', T cell receptor gene', 'immunoglobulin gene') 
#
#
############################################################################################################


def get_ppi(lcc):
    
    # Open database connection
    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabPW>","GenesGO")
    # prepare a cursor object using cursor() method
    cursor = db.cursor()

    sql = """
            SELECT
            e.entrez_1,
            e.entrez_2,
            g1.Locus_Type,
            g1.Locus_Group,
            g2.Locus_Type,
            g2.Locus_Group
            FROM networks.PPI_hippie2017 e
            INNER JOIN GenesGO.hgnc_complete g1 ON e.entrez_1 = g1.Entrez_Gene_ID_NCBI
            INNER JOIN GenesGO.hgnc_complete g2 ON e.entrez_2 = g2.Entrez_Gene_ID_NCBI
            WHERE 
                (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'T cell receptor gene')	          # 0 links    
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'immunoglobulin gene')              # 4 links    
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'T cell receptor gene')             # 0 links
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'immunoglobulin gene') 	          # 0 links         
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'T cell receptor gene' AND g2.Locus_Type = 'gene with protein product')       # 17 links         
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'T cell receptor gene')       # 1 links         
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'immunoglobulin gene' AND g2.Locus_Type = 'gene with protein product')        # 115 links         
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'immunoglobulin gene')        # 295 links    
                OR (e.author != '' AND e.entrez_1 != e.entrez_2
                    AND g1.Locus_Type = 'gene with protein product' AND g2.Locus_Type = 'gene with protein product')  # 309602 links  

        
        """
    try: 
        # execute SQL query 
        cursor.execute(sql)
        data = cursor.fetchall()

    except:
        print('SQL error')

    db.close()    
    
    l_nodes = []
    for x in data:
        l_nodes.append(x[0])
        l_nodes.append(x[1])
    l_nodes = list(set(l_nodes))

    G = nx.Graph()
    G.add_nodes_from(l_nodes)
    
    for x in data:
        G.add_edge(x[0],x[1])
    
    if lcc == 1:
        Nl_l = sorted(nx.connected_components(G)) # generates list of components node lists
        l = [len(x) for x in Nl_l] # returns list of length of node lists 
        idx = l.index(max(l))   # find the index of the maximal length i.e. lcc
        Nlcc = Nl_l[idx]    # pin down lcc
        G_lcc = G.subgraph(Nlcc)   # extract lcc graph
        G = G_lcc.copy()
    else:
        pass
    
    return G




######################
#GO

"""
function gets (UO or DOWN DIRECTED) GO branch and if it is (un-)directed
and returns nx.Graph
categories: Function, Process, Component

"""
def GO_tree(branch_category,graphtype,stream = 'up'):
    
    if branch_category == 'Component':
        cat = 'cellular_component'
    if branch_category == 'Function':
        cat = 'molecular_function'
    if branch_category == 'Process':
        cat = 'biological_process'        
    
    db = pymysql.connect("menchelabdb.int.cemm.at","readonly","ra4Roh7ohdee","GenesGO")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql_tree =  """
                    SELECT DISTINCT go_id,is_a,relationship 
                    FROM GO_tree WHERE is_a != '' 
                    AND namespace in ('%s')
                """ %(cat)
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql_tree)
        data_tree = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()
    if graphtype == 'directed':
        G = nx.DiGraph()
    else:
        G = nx.Graph()
        
    l_GOterms = []
    c = 0
    for x in data_tree:
        c +=1
        #if c < 10:
        #l_GOterms.append(x[0])
        G.add_node(x[0])

        
    if branch_category == 'Component':
        # 'cellular_component':'GO:0005575'
        tip = 'GO:0005575'
    if branch_category == 'Function':
        # 'molecular_function':'GO:0003674'  
        tip = 'GO:0003674'
    if branch_category == 'Process':
        # 'biological_process':'GO:0008150'  
        tip = 'GO:0008150'
    G.add_node(tip)

    if stream == 'up':
        
        # 1st condition to make a link
        # Node A get a directed link to B whenever there is B in the list of the 'is_a' column
        c = 0
        for x in data_tree:
            go_id = x[0]
            is_a = x[1]
            n_entries_inCol = int(is_a.count('!')) 
            for n in range(0,n_entries_inCol):
                c += 1
                G.add_edge(go_id,is_a.split('!')[n][-11:-1])

        # 2nd condition to make a link
        # Node A get a directed link to B whenever there is B as part_of/regulates in the list of the 'relationship' column
        c = 0
        for x in data_tree:
            go_id = x[0]
            is_a = x[1]
            rel = x[2]
            n_entries_inCol = int(rel.count('!')) 
            if n_entries_inCol > 0:
                for n in range(0,n_entries_inCol):
                    c += 1
                    #print(rel.split('!')[n][-11:-1])
                    G.add_edge(go_id,rel.split('!')[n][-11:-1])        

    else:
        # down   

        # 1st condition to make a link
        # Node A get a directed link to B whenever there is B in the list of the 'is_a' column
        c = 0
        for x in data_tree:
            go_id = x[0]
            is_a = x[1]
            n_entries_inCol = int(is_a.count('!')) 
            for n in range(0,n_entries_inCol):
                c += 1
                G.add_edge(is_a.split('!')[n][-11:-1],go_id)

        # 2nd condition to make a link
        # Node A get a directed link to B whenever there is B as part_of/regulates in the list of the 'relationship' column
        c = 0
        for x in data_tree:
            go_id = x[0]
            is_a = x[1]
            rel = x[2]
            n_entries_inCol = int(rel.count('!')) 
            if n_entries_inCol > 0:
                for n in range(0,n_entries_inCol):
                    c += 1
                    #print(rel.split('!')[n][-11:-1])
                    G.add_edge(rel.split('!')[n][-11:-1],go_id)
                
                

    return G, tip


"""

function gets list of GO term-IDs as input 
and returns a dictionary of keys: GO-IDs values: GO-names

"""

def turnGOids_GOnames(list_goterms):
    
    # Open database connection
    db = pymysql.connect("menchelabdb.int.cemm.at","readonly","ra4Roh7ohdee","GenesGO")
    # prepare a cursor object using cursor() method
    cursor = db.cursor()

    str_goID = ''
    c = 0
    # Read sample gene set here and make it SQL queryable 
    for gd in list_goterms:
        str_goID += '"' + str(gd) + '",'
        c += 1
    str_goID = str_goID[:len(str_goID)-1]  
    sql_goname = """ SELECT go_id,go_name FROM GO_tree WHERE go_id in (%s) """ %str_goID
    #cursor.execute(sql_goname)

    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql_goname)
        data_goidname = cursor.fetchall()

    except:
        print('load GO IDs to GO names -  SQL error')
    db.close()
    #print(data_goidname)
    d_go_id_name = {}    
    for x in data_goidname:
        d_go_id_name[x[0]] = (x[1])
    
    return d_go_id_name




"""
function gets GO to gene annotation from db
returns dictionary gene: go_list
and go: gene_list
categories: Function, Process, Component
"""
def loadgene2go(category = 'Function'):      
    
    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabPW>","GenesGO")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql =  """
                SELECT DISTINCT g.go_id, g.entrezid
                FROM GenesGO.Gene2GO_human g
                WHERE g.entrezid != '-'
                AND g.evidence != 'ND' 
                AND g.evidence != 'IEA' 
                AND g.evidence != 'IPI' 
                AND g.aspect  = '%s'
                """ %(category[0])
#     cursor.execute(sql)
#     data_gd = cursor.fetchall()
    
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql)
        data_gd = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()

    d_gene_go = defaultdict(list)
    d_go_gene = defaultdict(list)


    c = 0
    for xx in data_gd:
        c +=1
        gene = xx[1]
        go = xx[0]
        d_gene_go[gene].append(go)
        d_go_gene[go].append(gene)


    return d_gene_go, d_go_gene




"""
function returns dict key: go value: list_go
"""
def loadupstreams(category = 'Component'):      
    
    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabPW>","GenesGO")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql =  """
                SELECT
                go_id,
                up_ids
                FROM GenesGO.GO_upterms
                WHERE namespace = '%s'
                """ %(category[0])
#     cursor.execute(sql)
#     data_gd = cursor.fetchall()
    
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql)
        data_gd = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()

    d_go_up = defaultdict(list)

    c = 0
    for xx in data_gd:
        c +=1
        go = xx[0]
        up = xx[1]
        d_go_up[go].append(up)


    return d_go_up



####################################################################################

#Disease





"""
function gets DO tree and if it is (un-)directed
and returns nx.Graph

"""
def Disease_ontology(graphtype, stream = 'up'):      
    
    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabPW>","Gene2Disease")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql_tree =  """
                    SELECT 
                        do_id,
                        is_a-- ,
                    -- LENGTH(is_a) - LENGTH(REPLACE(is_a, '!', '')) counts
                    FROM Gene2Disease.disease_ontology
                    WHERE is_obsolete = ''
                    -- ORDER BY counts desc
                """ 
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql_tree)
        data_tree = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()
    if graphtype == 'directed':
        G = nx.DiGraph()
    else:
        G = nx.Graph()
        
    l_DOterms = []
    c = 0
    for x in data_tree:
        c +=1
        #if c < 10:
        #l_GOterms.append(x[0])
        G.add_node(x[0])

    # tip of dieseast tree: 'disease'
    tip = 'DOID:4'
    G.add_node('DOID:4')
    
    
    if stream == 'up':
        # 1st condition to make a link
        # Node A get a directed link to B whenever there is B in the list of the 'is_a' column
        c = 0
        for x in data_tree:
            do_id = x[0]
            is_a = x[1]
            n_entries_inCol = int(is_a.count('!')) 
#             print(n_entries_inCol)
            for n in range(1,n_entries_inCol+1):
                c += 1
                do_id_2 = str('DO' + str(is_a.split('DO')[n].split('!')[0])).strip()
#                 print(do_id)
#                 print(do_id_2[:-1], len(do_id_2))
                G.add_edge(do_id,do_id_2)
    else:
        # 1st condition to make a link
        # Node A get a directed link to B whenever there is B in the list of the 'is_a' column
        c = 0
        for x in data_tree:
            do_id = x[0]
            is_a = x[1]
            n_entries_inCol = int(is_a.count('!')) 
#             print(n_entries_inCol)
            for n in range(1,n_entries_inCol+1):
                c += 1
                do_id_2 = 'DO' + str(is_a.split('DO')[n].split('!')[0])
#                 print(do_id)
#                 print(do_id_2[:-1], len(do_id_2))
                G.add_edge(do_id_2[:-1],do_id)        
            

    return G,tip

"""
function gets DO to gene annotation from db
returns dictionary gene: do
"""
def diseaseID2name():      
    
    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabPW>","Gene2Disease")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql =  """
            SELECT
                do_id,
                do_name
            FROM disease_ontology
            """ 
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql)
        data = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()

    d_ID_name = {}
    for x in data:
        doID = x[0]
        do_name = x[1]
        d_ID_name[doID] = do_name


    return d_ID_name



"""
function gets DO to gene annotation from db
returns dictionary gene: do
"""
def loadgene2do():      

    #new version is called gene2disease_all

    db = pymysql.connect("<MencheLabServer>","readonly","<MencheLabServer>","Gene2Disease")
    cursor = db.cursor()
    # Query to contruct GO Hierarchy
    #
    #          Returns complete list of GO IDs that have a 'is a' relation 
    #
    sql =  """
                SELECT
                    mp.vocabulary,
                    mp.identifier,
                    gd.entrezId
                FROM Gene2Disease.gene2disease gd
                INNER JOIN HumanPhenotypes.umlsmapping mp
                ON mp.disease_ID = gd.diseaseID
                WHERE mp.vocabulary = 'DO'
                """ 
    try: 
        # execute SQL query using execute() method.
        cursor.execute(sql)
        data_gd = cursor.fetchall()
    except:
        print('GO tree SQL error')
    db.close()

    d_gene_do = defaultdict(list)
    d_do_gene = defaultdict(list)


    c = 0
    for xx in data_gd:
        c +=1
        gene = xx[2]
        do = str(xx[0]) + 'ID:' + str(xx[1]) 
        d_gene_do[gene].append(do)
        d_do_gene[do].append(gene)


    return d_gene_do, d_do_gene





def getAllGene_Annotation(go_branch):
    # CPU TIME BENCHMARK: FOR COMPLETE PPI AND GO-BRANCH BIOLOGICAL PROCESS RUNTIME APPROX 30 SEC

    d_gene_go, d_go_gene = loadgene2go(category = go_branch)

    d_go_up = loadupstreams(go_branch)


    GO_genes_annotation = {}

    # GENERATE DICT - KEYS: GENES VALUES: LIST OF GO TERMS INCLUDING UPSTREAM
    d_gene_go_up = defaultdict(list)

    for gene,l_t in d_gene_go.items():

        l_d = []
        for t in l_t:
            l_d.append(d_go_up[t])
            l_d.append([t])


        l_d_flat = set([item for sublist in l_d for item in sublist])
        d_gene_go_up[gene] = l_d_flat




        for d in l_d_flat:
            if GO_genes_annotation.has_key(d):
                GO_genes_annotation[d].append(gene)
            else:
                GO_genes_annotation[d] = [gene]

    return d_gene_go,GO_genes_annotation


def getAllGene_Disease_Annotation():


    # LOAD GENE 2 DISEASE and vice versa
    d_gene_do, d_do_gene = loadgene2do()

    # LOAD ONTOLOGY TREE  HERE
    G_tree, tip = Disease_ontology( 'directed', 'up')

    # GENERATE DICT - KEYS: GENES VALUES: LIST OF DISEASE TERMS TERMS INCLUDING UPSTREAM
    d_nodes_terms_up = defaultdict(list)
    d_diseases_annotation = {}


    #names = diseaseID2name()

    #n = gene
    #l_t = leave disease terms
    l_exc = []
    for n, l_t in d_gene_do.items():


        # l_d = upstream genes
        l_d = []
        for t in l_t:
            try:
                for updos in nx.shortest_simple_paths(G_tree, t, tip):
                    for d in updos:
                        l_d.append(d)
            except:
                l_exc.append(t)

        for d in set(l_d):
            if d_diseases_annotation.has_key(d):
                d_diseases_annotation[d].append(n)
            else:
                d_diseases_annotation[d] = [n]

        d_nodes_terms_up[n] = set(l_d)

    return d_nodes_terms_up,d_diseases_annotation

