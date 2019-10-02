import numpy as np
from os import listdir
from os.path import isfile, join
import os
from matplotlib import  pylab as plt
from scipy import  special

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


# Define Experiment
table = 'DPN1018Batch2Per_Image'


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



# Effect size
def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(
        ((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)


def calculate_ZFactor(drug, dmso):
    '''
    print '---'
    print np.std(drug)
    print np.std(dmso)
    print np.mean(drug)
    print np.mean(dmso)
    print '---'
    '''

    return 1 - ((3 * np.std(drug) + 3 * np.std(dmso)) / (abs(np.mean(drug) - np.mean(dmso))))


def calculate_ZScore(drug_mean, dmso):
    '''

    '''

    return (drug_mean - np.mean(dmso))/np.std(dmso)


def calculate_Effect_size(table):
    '''
    This function calculates the effect size of a given feature using cohen's calculation. The higher the

    :return:
    '''

    ensure_dir('../results/'+table+'/Effect_Size/Effect_Sizes.csv')
    fp_out = open('../results/'+table+'/Effect_Size/Effect_Sizes.csv','w')
    fp_out.write("Feature,Single_Cohen's_D_(MedianTop5%),Single_ZScore_(MedianTop5%),Single_Number_Significant,Combination_ZScore_(MedianTop5%),Combination_Number_Significant\n")


    features = getFeatureList()
    single_significant = []
    comb_significant = []
    for f in features:
        print f

        screen_results = get_feature_result(f, table)
        plates = screen_results.keys()

        dmso_poc = []
        drug_poc = {}
        combination_poc = {}

        for p in plates:
            for well in screen_results[p]:
                if screen_results[p][well]['Worked'] == 'FALSE' or screen_results[p][well]['N_Value'] == -100000:
                    continue

                if screen_results[p][well]['Drug_1'] == 'DMSO':
                    dmso_poc.append(screen_results[p][well]['N_Value'])
                elif screen_results[p][well]['Drug_2'] == 'DMSO':
                    if drug_poc.has_key(screen_results[p][well]['Drug_1']):
                        drug_poc[screen_results[p][well]['Drug_1']].append(screen_results[p][well]['N_Value'])
                    else:
                        drug_poc[screen_results[p][well]['Drug_1']] = [screen_results[p][well]['N_Value']]
                elif 'CLOUD' in screen_results[p][well]['Drug_2']:
                    combination_poc[screen_results[p][well]['Drug_1']+','+screen_results[p][well]['Drug_2']] = screen_results[p][well]['N_Value']

        dmso_poc = reject_outliers_2(dmso_poc)  # DMSO values
        single_drug_cohenDs = []
        single_drug_ZScore = []
        combination_ZScore = []
        for key in drug_poc:
            drug_vals = reject_outliers_2(drug_poc[key])  # drug values
            if len(drug_poc[key]) > 1:
                # calculate cohens D between single drug values and DMSO mean
                cd = abs(cohen_d(drug_vals, dmso_poc))  # abs = absolute (turns negative values into positive)
                single_drug_cohenDs.append(cd)
            if len(drug_poc[key]) >= 1:
                single_drug_ZScore.append(abs(calculate_ZScore(np.mean(drug_vals),dmso_poc)))
            else:
                continue

        for comb_val in combination_poc.values():

            ZScore = abs(calculate_ZScore(comb_val,dmso_poc))
            combination_ZScore.append(ZScore)

        p_values_combination_ZScore  = (1 -  special.ndtr(combination_ZScore))*len(combination_ZScore)
        for (i, item) in enumerate(p_values_combination_ZScore):
            if item > 1:
                p_values_combination_ZScore[i] = 1

        p_values_single_drug_ZScore = (1 - special.ndtr(single_drug_ZScore)) * len(single_drug_ZScore)
        for (i, item) in enumerate(p_values_single_drug_ZScore):
            if item > 1:
                p_values_single_drug_ZScore[i] = 1



        fp_out.write(f+','+str(np.median(np.percentile(single_drug_cohenDs,95)))+
                     ','+str(np.median(np.percentile(single_drug_ZScore,95)))+
                     ','+str(len([x for x in p_values_single_drug_ZScore if x < 0.05])) +
                     ',' + str(np.median(np.percentile(combination_ZScore, 95))) +
                     ',' + str(len([x for x in p_values_combination_ZScore if x < 0.05])) + '\n')

        single_significant.append(len([x for x in p_values_single_drug_ZScore if x < 0.05]))
        comb_significant.append(len([x for x in p_values_combination_ZScore if x < 0.05]))

    fp_out.close()


    plt.hist(single_significant,bins='auto',color='grey')
    plt.axvline(1,ls='--',c='red')
    plt.legend(['Rejected Features: %d' %len([x for x in single_significant if x < 1])] )
    plt.savefig('../results/'+table+'/Effect_Size/Single_Effect_Sizes.pdf')


    plt.hist(comb_significant,bins='auto',color='grey')
    plt.axvline(1,ls='--',c='red')
    plt.legend(['Rejected Features: %d' %len([x for x in comb_significant if x < 1])] )
    plt.savefig('../results/'+table+'/Effect_Size/Combination_Effect_Sizes.pdf')


print 'Calculate Effect size for: %s' %table

calculate_Effect_size(table)




