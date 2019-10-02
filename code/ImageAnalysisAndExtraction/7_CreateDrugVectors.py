import numpy as np
import pickle
from matplotlib import pylab as plt
import os
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA
import matplotlib
from scipy.spatial import distance
from os import listdir
from os.path import isfile, join
import numpy as np
import umap
import random
from sklearn.manifold import TSNE


# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
table = 'DPN1018Batch1Per_Image'
table2 = 'DPN1018Batch2Per_Image'


# - - - - - - - - - - - - - - - - - - - -


def Plot_PCA(x_val,col,y_val='None'):
    '''
    create pca plot

    :param x_val: n-dimensional feature string for each drug
    :param y_val: drug name
    :return:
    '''
    drug_x = x_val
    if y_val != "None":
        drug_name=y_val




    pca = PCA(n_components=2)
    pca.fit(drug_x)
    print pca.explained_variance_

    X = pca.transform(drug_x)
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], alpha=0.4,c=col)

    if y_val != "None":
        for i, txt in enumerate(drug_name):
            ax.annotate(txt, (X[:, 0][i], X[:, 1][i]), size=8)

    # ax.set_xlabel(str(pca.explained_variance_[0]))
    # ax.set_ylabel(str(pca.explained_variance_[1]))

    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')

    plt.show()
    #plt.savefig('../results/Batch' + str(batch) + '/Remove_Replicates/PCA.pdf', format='pdf', dpi=600)
    plt.close()

#Some Easy Outlier detection
def reject_outliers_2(data, m = 6.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return [data[i] for i in range(0,len(data)) if s[i] < m]

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




# ---------


def make_DrugVectors(table,UseCombinedFeatures):
    '''
    Function to create the unified, as well as individual single treatment vectors


    :param batch:
    :return:
    '''


    #Get allowed features
    features = []

    if UseCombinedFeatures:
        fp = open('../results/' + table + '/Remove_Correlation/Uncorrelated_Features_CombinedBatches.csv', 'r')
    else:
        fp = open('../results/' + table + '/Remove_Correlation/Uncorrelated_Features.csv', 'r')

    for line in fp:
        features.append(line.strip())
    print 'Number of ok Features: %d' % len(features)
    fp.close()
    features.sort()
    features = features


    treatments = {}
    for f in features:
        print f
        screen_results = get_feature_result(f,table)
        plates = screen_results.keys()



        #Find for each plate the corresponding SINGLE, and calculate a feature vector as the difference between
        #DMSO and the actual value
        for p in plates:
            dmso_poc = []
            all_Wells = {}

            for well in screen_results[p]:



                if screen_results[p][well]['Worked'] == 'FALSE' or screen_results[p][well]['N_Value'] == -100000:
                    continue

                if screen_results[p][well]['Drug_1'] == 'DMSO':
                    dmso_poc.append(screen_results[p][well]['N_Value'])




                #drug_poc[screen_results[p][well]['Drug_1']+','+screen_results[p][well]['Drug_2']] = screen_results[p][well]['N_Value']
                all_Wells[screen_results[p][well]['Drug_1']+'|'+screen_results[p][well]['Drug_2']+'_'+well+'_'+str(p)] = screen_results[p][well]['N_Value']


            dmso_poc = reject_outliers_2(dmso_poc)
            mean_dmso_effect = np.mean(dmso_poc)


            for key in all_Wells:
                if treatments.has_key(key):
                    treatments[key].append(all_Wells[key] - mean_dmso_effect)
                else:
                    treatments[key] = [all_Wells[key] - mean_dmso_effect]


    all_treatments = treatments.keys()
    all_treatments.sort()

    ensure_dir('../results/'+table+'/PerturbationVectors/Vectors.csv')

    if UseCombinedFeatures:
        fp_out = open('../results/' + table + '/PerturbationVectors/Vectors_Combined.csv', 'w')
    else:
        fp_out = open('../results/'+table+'/PerturbationVectors/Vectors.csv', 'w')
    fp_out.write('Perturbation,' + ','.join(features) + '\n')



    for t in all_treatments:
        fp_out.write(t+','+','.join([str(treatments[t][x]) for x in range(0,len(features))])+'\n')
    fp_out.close()


make_DrugVectors(table2,True)
#make_DrugVectors(table2,True)



def make_combined_VectorList(table1,table2):



    if os.path.isfile('../results/' + table1 + '/PerturbationVectors/Vectors_Combined.csv') == False or os.path.isfile('../results/' + table2 + '/PerturbationVectors/Vectors_Combined.csv') == False:
        print 'First create single files'
        exit()
    else:



        fp_out = open('../results/All_Vectors_Combined.csv','w')


        fp = open('../results/' + table1 + '/PerturbationVectors/Vectors_Combined.csv')
        for line in fp:
            fp_out.write(line)
        fp.close()


        fp = open('../results/' + table2 + '/PerturbationVectors/Vectors_Combined.csv')
        fp.next()
        for line in fp:
            fp_out.write(line)
        fp.close()

        fp_out.close()

#make_combined_VectorList(table,table2)




