import networkx as nx
import numpy as np
from os import listdir
from os.path import isfile, join
from matplotlib import pylab as plt
from collections import Counter
import os
import seaborn as sns


# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
table = 'DPN1018Batch1Per_Image'
table2 = 'DPN1018Batch2Per_Image'
# - - - - - - - - - - - - - - - - - - - -

# Some Easy Outlier detection
def reject_outliers_2(data, m=6.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / (mdev if mdev else 1.)
    return [data[i] for i in range(0, len(data)) if s[i] < m]

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)



# Methods for getting feature by feature
def getFeatureList(mypath='../results/' + table + '/POCNormalized/'):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyfiles.sort()
    features = []
    for file in onlyfiles:
        features.append(file.strip().split('.')[0])


    if '' in features:
        features.remove('')

    if 'MaxMin_Values' in features:
        features.remove('MaxMin_Values')


    return features

def get_feature_result(feature, db_table):

    # go through the input file of the feature
    # path = '../data/Normalized_Wells/' + feature + '.csv'
    path = '../results/' + db_table + '/POCNormalized/' + feature + '.csv'
    fp = open(path, 'r')
    fp.next()
    feature_results = {}
    # mean = {}
    for line in fp:
        tmp = line.strip().split(',')

        plate = int(tmp[0])
        well = tmp[1]
        drug1 = tmp[2]
        drug2 = tmp[3]
        worked = tmp[4]

        # if 'nan' for some features this might happen, then just set to mean of the plate
        if tmp[5] != 'nan':
            normed_value = tmp[5]

        else:
            # normed_value = np.mean(mean[plate])
            normed_value = 0
            worked = 'FALSE'

        if normed_value == -100000.0:
            worked = 'FALSE'


        #else basically create an entry in the dictionary with the information as well as the normed value
        #if the dictionary does not yet contain the plate, then set it
        if feature_results.has_key(plate):
            feature_results[plate][well] = {'Drug_1': drug1, 'Drug_2': drug2, 'Worked': worked,
                                                 'N_Value': float(normed_value)}
        else:
            feature_results[plate] = {
                well: {'Drug_1': drug1, 'Drug_2': drug2, 'Worked': worked, 'N_Value': float(normed_value)}}

    # print feature_results
    return feature_results

def get_FilteredFeatures(table):
    '''
    Filter features

    :param table:
    :return:
    '''

    ###
    # IntraCV
    ###
    fp = open('../results/' + table + '/Coefficient_of_Variation/IntraPlate_Variability.csv', 'r')
    fp.next()
    intra_CV_features_passed = set()
    for line in fp:
        tmp = line.strip().split(',')
        if float(tmp[1]) < 0.2:
            intra_CV_features_passed.add(tmp[0])
    fp.close()

    ###
    # InterCV
    ###
    fp = open('../results/' + table + '/Coefficient_of_Variation/InterPlate_Variability.csv', 'r')
    fp.next()
    inter_CV_features_passed = set()
    for line in fp:
        tmp = line.strip().split(',')
        if float(tmp[1]) < 0.2:
            inter_CV_features_passed.add(tmp[0])
    fp.close()

    ###
    # Replicate Correlation (Reproduceability)
    ###
    fp = open('../results/' + table + '/Correlation_Between_Replicates/Correlation_Results.csv', 'r')
    fp.next()
    replicate_Cor_features_passed = set()
    for line in fp:
        tmp = line.strip().split(',')
        if float(tmp[1]) > 0.2:
            replicate_Cor_features_passed.add(tmp[0])
    fp.close()


    ###
    # Effect Size
    ###
    fp = open('../results/' + table + '/Effect_Size/Effect_Sizes.csv', 'r')
    fp.next()
    effect_size_features_passed = set()
    for line in fp:
        tmp = line.strip().split(',')
        if float(tmp[5]) > 0:
            effect_size_features_passed.add(tmp[0])
    fp.close()


    all_passed = intra_CV_features_passed.intersection(inter_CV_features_passed).intersection(replicate_Cor_features_passed).intersection(effect_size_features_passed)

    print 'Features passed (%s):' %table
    print 'IntraCV: %d' % len(intra_CV_features_passed)
    print 'InterCV: %d' %len(inter_CV_features_passed)
    print 'Replicate Correlation: %d' % len(replicate_Cor_features_passed)
    print 'Effect Size: %d' % len(effect_size_features_passed)
    print 'Union: %d' %len(all_passed)
    print '----'
    return all_passed


def create_correlation_Network_AllFeatures(table, saveNetwork=False):


    fp = open('../results/' + table + '/Effect_Size/Effect_Sizes.csv', 'r')
    fp.next()
    effect_size_dic = {}
    for line in fp:
        tmp = line.strip().split(',')

        effect_size_dic[tmp[0]] = int(tmp[5])
    fp.close()


    features = getFeatureList()


    feature_dmso_values = {}

    for f in features:
        screen_results = get_feature_result(f, table)
        plates = screen_results.keys()

        dmso_values = []
        for plate in plates:
            for well in screen_results[plate]:
                if screen_results[plate][well]['Worked'] == 'FALSE' or screen_results[plate][well][
                    'N_Value'] == -100000:
                    continue
                if screen_results[plate][well]['Drug_1'] == 'DMSO':
                    dmso_values.append(screen_results[plate][well]['N_Value'])
        feature_dmso_values[f] = dmso_values

    # create output graph
    G = nx.Graph()

    # 1.) go through all pairwise correlations
    # 2.) Calculate pearson correlation
    # 3.) If correlation bigger than 0.7 create edge in output graph

    # print 'Calculate Correlation'
    heatmap_data = []
    for f in feature_dmso_values.keys():
        print f
        feature1 = feature_dmso_values[f]
        G.add_node(f)
        tmp = []
        for f2 in feature_dmso_values.keys():
            #if f != f2:
            feature2 = feature_dmso_values[f2]
            cor = np.corrcoef(feature1, feature2)[0, 1]

            if abs(cor) > 1:
                G.add_edge(f, f2, weight=cor)
                G.node[f]['Effect'] = effect_size_dic[f]
                G.node[f2]['Effect'] = effect_size_dic[f2]
            tmp.append(cor)
        heatmap_data.append(tmp)

    # Save correlation network for manual inspection
    if saveNetwork:
        ensure_dir('../results/'+table+'/Remove_Correlation/CorrelationNetwork.gml')
        nx.write_gml(G, '../results/'+table+'/Remove_Correlation/CorrelationNetwork.gml')
        sns.clustermap(data=heatmap_data, cmap="RdBu")
        #plt.show()
        plt.savefig('../results/'+table+'/Remove_Correlation/CorrelationHeatMap.pdf')
        plt.close()

    return G

# Correlation Functions
def check_edge_between(min_nod, g):
    '''
    Find minimal amount of nodes to remove, so that min_nodes are not neighbors anymore

    :param min_nod: list of nodes to check
    :param g: connected component subgraph
    :return: list of nodes to remove so that nodes in min_nod are not connected
    '''

    # extract all pairwise edges between nodes in min_nod
    has_edge = []
    for min1 in min_nod:
        for min2 in min_nod:
            if min1 > min2:
                if g.has_edge(min1, min2):
                    has_edge.append(min1)
                    has_edge.append(min2)

    # create Counter instance
    data = Counter(has_edge)
    # get values of which node occured how often edges
    # e.g. [2,2,1,1,1,1] = two nodes are connected with two in min_nodes while 4 nodes are connected only with one
    freq_list = data.values()

    # if freq_list == 0, all nodes are separated from each other
    if len(freq_list) == 0:
        return []

    # find max edges of on of the nodes (e.g. example above would be 2)
    max_cnt = max(freq_list)

    # get all nodes that are involved in max_cnt (e.g. the first two nodes from the freq_list)
    total = freq_list.count(max_cnt)

    # if all nodes are equaly e.g. [1,1,1], it does'nt bother which one to remove, choose randomly one (so remove the other two)
    if total == len(freq_list):
        max_val = 0
        max_node = ''
        for node in min_nod:
            if g.node[node]['Effect'] > max_val:
                max_val = g.node[node]['Effect']
                max_node = node

        keep = max_node
        # keep = choice(min_nod)
        copylist = list(min_nod)
        copylist.remove(keep)
        return copylist

    # Return these nodes for removal (first two from example above)
    most_common = data.most_common(total)
    return [elem[0] for elem in most_common]


def analyse_component(g, draw=False):
    '''
    Main function for max fragmentation of subgraphs.
    Takes a graph (origins from bigger network as connected component)
    Slowly fragmentises it by removing best suited nodes

    :param g: connected component subgraph
    :param draw: True if there should be a step by step drawing output
    :return: number of
    '''
    # contains max number nodes for current component
    tmp_keep = []

    # fragmentise connected component subgraph until all nodes fragmented
    while len(g.nodes()) > 0:

        # draw option
        if draw == True:
            nx.draw_networkx(g, pos=nx.spring_layout(g), with_labels=False)
            plt.draw()  # pyplot draw()
            plt.show()

        # list for nodes that need to be removed in each iteration
        # Contains: Selected Nodes (find in tmp_keep), as well as their neighbors
        nodes_to_remove = set()

        # if (remaing) component is only two nodes ==>  A--B; take randomly one of the two
        if len(g.nodes()) == 2 and len(g.edges()) == 1:
            two_nodes = list(g.nodes())

            if g.node[two_nodes[0]]['Effect'] > g.node[two_nodes[1]]['Effect']:
                tmp_keep.append(two_nodes[0])
            else:
                tmp_keep.append(two_nodes[1])

            # purely random choice of which node to keep
            # rand_node = choice(g.nodes())
            # tmp_keep.append(rand_node)

            nodes_to_remove.add(list(g.nodes())[0])
            nodes_to_remove.add(list(g.nodes())[1])

        # if bigger than only two connected nodes

        else:
            # get node degrees
            degrees_tmp = g.degree()

            degrees = {}
            for d in degrees_tmp:
                degrees[d[0]] = d[1]

            # find terminal nodes (= degree 1)
            terminal_nodes = [x for x in degrees if degrees[x] == 1]

            # if subgraph still has terminal nodes, choose these
            if len(terminal_nodes) > 0:
                for tn in terminal_nodes:
                    tmp_keep.append(tn)
                    nodes_to_remove.add(list(g.edges(tn))[0][1])
                    nodes_to_remove.add(tn)

            # if no terminal nodes exist
            else:

                # Check if there are nodes with higher degree than other
                # if all degrees uniformly it's for example a triangle, rectangle etc. (circularity)
                if all(x == degrees.values()[0] for x in degrees.values()) == False:
                    # example for some nodes with lower/higher degree than others
                    # A-B
                    # |\|  ==> in this case the algorithm should pick B and C (other alternative would be only A or D)
                    # C-D

                    # extract smalles degree
                    min_degree = min(degrees.values())

                    # get nodes with this smallest degree
                    min_nodes = [x for x in degrees if degrees[x] == min_degree]

                    # check if these nodes with smallest degree are somehow neighbors
                    # e.g. two rectangles (4 nodes) connected by a middle node
                    # ==> always the three "outer" rectangle nodes would have degree 2 (togher with the middle one connecting the two rectangles)
                    while True:
                        # remove the minimum amount of nodes, so all selected "min_nodes" are no neighbors anymore
                        node_edge_remove = check_edge_between(min_nodes, g)
                        if len(node_edge_remove) == 0:
                            break
                        for node in node_edge_remove:
                            min_nodes.remove(node)

                    # Save the Min_nodes to tmp_keep and add them to nodes_to_remove (togher with their neighbors)
                    for mn in min_nodes:
                        tmp_keep.append(mn)
                        nodes_to_remove.add(mn)
                        edges = g.edges(mn)
                        for edge in edges:
                            nodes_to_remove.add(edge[1])

                # if all degrees are uniformly, meaning you have a triangle, rectangle, fully connected graph
                # e.g.:
                # A-B
                # | |  ==> e.g. first pick A (randomly); remove B + C (neighbors); in next
                # C-D           iteration there is only D left (will be picked)
                else:
                    # randomly choose a single node (all nodes equally anyway)
                    max_val = 0
                    max_node = ''
                    for node in g.nodes():
                        if g.node[node]['Effect'] > max_val:
                            max_val = g.node[node]['Effect']
                            max_node = node

                    rand_node = max_node

                    # rand_node = choice(g.nodes())
                    # add this random nood to tmp_keep and again remove him together with the neighbors
                    tmp_keep.append(rand_node)
                    nodes_to_remove.add(rand_node)
                    edges = g.edges(rand_node)
                    for edge in edges:
                        nodes_to_remove.add(edge[1])

        # Remove nodes from current subgraph
        for ntr in nodes_to_remove:
            g.remove_node(ntr)

        if draw:
            print tmp_keep

    return tmp_keep


def remove_Correlating_Features(table,corelation_Threshold=0.6):

    filtered_features = get_FilteredFeatures(table)


    print 'Create Correlation network for: %s' %table
    # get network file
    if os.path.isfile('../results/'+table+'/Remove_Correlation/CorrelationNetwork.gml') == True:
        print 'Using already existing correlation network!'
        g = nx.read_gml('../results/'+table+'/Remove_Correlation/CorrelationNetwork.gml')
    else:
        g = create_correlation_Network_AllFeatures(table, saveNetwork=True)

    exit()

    #Remove self edges (by definitation cor = 1)
    g.remove_edges_from(g.selfloop_edges())

    # Remove non correlating edges
    all_edges =  list(g.edges())
    for edge in all_edges:
        cor = abs(g[edge[0]][edge[1]]['weight'])
        if cor < corelation_Threshold:
            g.remove_edge(edge[0],edge[1])

    all_nodes = list(g.nodes())
    for node in all_nodes:
        if node not in filtered_features:
            g.remove_node(node)


    print 'Correlation Network (All Features): '
    print 'Number of nodes: %d' %len(g.nodes())
    print 'Number of edges: %d' %len(g.edges())
    print '--'





    # Extract connected components
    components = nx.connected_component_subgraphs(g)

    # keep contains the max amount of nodes that are never connected
    keep = []

    # Go threw components
    for comp in components:
        # if single node, add to keep (anyway not connected to anything)
        if len(comp.nodes()) == 1:
            keep.append(list(comp.nodes())[0])
        # if not single node check maximum fragmentation
        else:
            keep = keep + analyse_component(comp, False)



    keep.sort()

    ensure_dir('../results/'+table+'/Remove_Correlation/Uncorrelated_Features.csv')
    fp = open('../results/'+table+'/Remove_Correlation/Uncorrelated_Features.csv','w')

    for feature in keep:
        fp.write(feature+'\n')
    fp.close()
    print 'Number of Features: %d' %len(keep)


#remove_Correlating_Features(table2)



def remove_Correlating_Features_CombineBatches(table,table2,corelation_Threshold=0.6):



    filtered_features1 = get_FilteredFeatures(table)
    filtered_features2 = get_FilteredFeatures(table2)

    if os.path.isfile('../results/' + table + '/Remove_Correlation/CorrelationNetwork.gml') == True and os.path.isfile('../results/' + table2 + '/Remove_Correlation/CorrelationNetwork.gml') == True:
            g1 = nx.read_gml('../results/' + table + '/Remove_Correlation/CorrelationNetwork.gml')
            g2 = nx.read_gml('../results/' + table2 + '/Remove_Correlation/CorrelationNetwork.gml')
    else:
        print 'First analyse batches separately'
        exit()




    #Remove self edges (by definitation cor = 1)
    g1.remove_edges_from(g1.selfloop_edges())
    g2.remove_edges_from(g2.selfloop_edges())

    all_nodes = list(g1.nodes())
    for node in all_nodes:
        if node not in filtered_features1 or node not in filtered_features2:
            g1.remove_node(node)


    # Remove non correlating edges
    all_edges =  list(g1.edges())
    for edge in all_edges:
        cor1 = abs(g1[edge[0]][edge[1]]['weight'])
        cor2 = abs(g2[edge[0]][edge[1]]['weight'])

        if cor1 < corelation_Threshold and cor2 < corelation_Threshold:
            g1.remove_edge(edge[0],edge[1])


    print 'Correlation Network (All Features): '
    print 'Number of nodes: %d' %len(g1.nodes())
    print 'Number of edges: %d' %len(g1.edges())
    print '--'


    # Extract connected components
    components = nx.connected_component_subgraphs(g1)

    # keep contains the max amount of nodes that are never connected
    keep = []

    # Go threw components
    for comp in components:
        # if single node, add to keep (anyway not connected to anything)
        if len(comp.nodes()) == 1:
            keep.append(list(comp.nodes())[0])
        # if not single node check maximum fragmentation
        else:
            keep = keep + analyse_component(comp, True)



    keep.sort()

    ensure_dir('../results/'+table+'/Remove_Correlation/Uncorrelated_Features_CombinedBatches.csv')
    ensure_dir('../results/' + table2 + '/Remove_Correlation/Uncorrelated_Features_CombinedBatches.csv')

    fp = open('../results/'+table+'/Remove_Correlation/Uncorrelated_Features_CombinedBatches.csv','w')
    fp2 = open('../results/' + table2 + '/Remove_Correlation/Uncorrelated_Features_CombinedBatches.csv', 'w')

    for feature in keep:
        fp.write(feature+'\n')
        fp2.write(feature + '\n')
    fp.close()
    fp2.close()
    print 'Number of Features: %d' %len(keep)


#remove_Correlating_Features_CombineBatches(table,table2,0.8)