import numpy as np
from os import listdir
from os.path import isfile, join
import os
from matplotlib import pylab as plt


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


# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
table = 'DPN1018Batch2Per_Image'


# - - - - - - - - - - - - - - - - - - - -


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


def calculate_intraPlate_CV(table):
    '''
    Calculate intra plate cv, should be smaller than 0.1. Bigger values indicate features that fluctuate too much
    in the same plate. Create a IntraPlate_CV file
    :param batch: which batch should be analysed
    :return: nothing
    '''
    print 'Calculate Intra Plate CV for plates (based on DMSO wells): %s' % (str(table))

    # load all features
    features = getFeatureList()

    # open output file
    ensure_dir('../results/'+table+'/Coefficient_of_Variation/IntraPlate_Variability.csv')
    fp_out = open('../results/'+table+'/Coefficient_of_Variation/IntraPlate_Variability.csv', 'w')
    fp_out.write('Feature,IntraPlateCV\n')


    intra_CVs = []
    # Calculate CV for each feature
    for f in features:
        # print f
        # load the data for the specific feature
        screen_results = get_feature_result(f, table)
        plates = screen_results.keys()
        # plate_cvs contains the individual intraplate cv's, (final cv = mean over all)
        plate_cvs = []

        for plate in plates:
            # plate well contains all the dmso (POC) values
            plate_wells = []

            # go through all the wells and look for DMSO treated wells
            for well in screen_results[plate]:
                if screen_results[plate][well]['Worked'] == 'FALSE' or screen_results[plate][well][
                    'N_Value'] == -100000:
                    continue

                if screen_results[plate][well]['Drug_1'] == 'DMSO' and screen_results[plate][well]['Drug_2'] == 'None':
                    plate_wells.append(screen_results[plate][well]['N_Value'])

            # if all valus below zero (now set to zero, means all dmso again same value) good cv, else divide by zero
            if sum(plate_wells) == 0:
                plate_cvs.append(0)
            else:
                # Reject most extreme outlier
                plate_wells_no_outlier = reject_outliers_2(plate_wells)

                # calculate the cv
                sd = abs(np.std(plate_wells_no_outlier))
                mean = abs(np.mean(plate_wells_no_outlier))
                cv = (sd / mean)

                # add the cv of this plate to the other plates
                plate_cvs.append(cv)

        # calculate mean cv over all plates and write to file
        fp_out.write(f + ',' + str(np.mean(plate_cvs)) + '\n')
        intra_CVs.append(np.mean(plate_cvs))

    fp_out.close()

    plt.hist(intra_CVs,bins='auto',color='grey')
    plt.axvline(0.2,ls='--',c='red')
    plt.legend(['Rejected Features: %d' %len([x for x in intra_CVs if x > 0.2])] )
    plt.savefig('../results/'+table+'/Coefficient_of_Variation/IntraPlate_Variability.pdf')
    plt.close()

#calculate_intraPlate_CV(table)


def calculate_interPlate_CV(table):
    '''
    Calculate intra plate cv, should be smaller than 0.1. Bigger values indicate features that fluctuate too much
    in the same plate. Create a IntraPlate_CV file
    :param batch: which batch should be analysed
    :return: nothing
    '''
    print 'Calculate Inter Plate CV for plates (based on DMSO wells): %s' % (str(table))

    # load all features
    features = getFeatureList()

    # open output file
    ensure_dir('../results/' + table + '/Coefficient_of_Variation/InterPlate_Variability.csv')
    fp_out = open('../results/' + table + '/Coefficient_of_Variation/InterPlate_Variability.csv', 'w')

    fp_out.write('Feature,IntraPlateCV\n')

    inter_CVs = []
    # Calculate CV for each feature
    for f in features:
        # load the data for the spcific feature
        screen_results = get_feature_result(f, table)
        plates = screen_results.keys()

        # plate_means means of all DMSO wells on one plate
        plate_means = []
        for plate in plates:
            # plate well contains all the dmso (POC) values
            plate_wells = []
            # go threw all the wells and look for DMSO treated wlls
            for well in screen_results[plate]:
                if screen_results[plate][well]['Worked'] == 'FALSE' or screen_results[plate][well][
                    'N_Value'] == -100000:
                    continue

                if screen_results[plate][well]['Drug_1'] == 'DMSO' and screen_results[plate][well][
                        'Drug_2'] == 'None':
                    plate_wells.append(screen_results[plate][well]['N_Value'])

            # Reject most extreme outlier
            plate_wells_no_outlier = reject_outliers_2(plate_wells)
            # add the mean of this plate, later calculate cv for the individual plates

            plate_means.append(np.mean(plate_wells_no_outlier))


        # if all values are zero, that means also pretty consistant (but would lead to divide by zero error)
        if sum(plate_means) == 0:
            cv = 0
        else:
            # calculate the cv
            sd = abs(np.std(plate_means))
            mean = abs(np.mean(plate_means))
            cv = (sd / mean)

        # write inter plate CV to file
        fp_out.write(f + ',' + str(cv) + '\n')
        inter_CVs.append(cv)


    fp_out.close()


    plt.hist(inter_CVs, bins='auto', color='grey')
    plt.axvline(0.2, ls='--', c='red')
    plt.legend(['Rejected Features: %d' % len([x for x in inter_CVs if x > 0.2])])
    plt.savefig('../results/' + table + '/Coefficient_of_Variation/InterPlate_Variability.pdf')
    plt.close()
    # print 'DONE'
#calculate_interPlate_CV(table)


def calculate_Correlation_BetweenPlates(table,plot=False):


    print 'Calculate Correlation between replicates (based on single drug wells): %s' % (str(table))

    # load all features
    features = getFeatureList()

    ensure_dir('../results/' + table + '/Correlation_Between_Replicates/')
    fp = open('../results/' + table + '/Correlation_Between_Replicates/Correlation_Results.csv','w')
    fp.write('Feature,PearsonCorrelation\n')
    correlations = []
    for f in features:

        #ensure_dir('../results/' + table + '/Replicates/' + f + '.csv')
        # fp_replicates = open('../results/Replicates_Batch' + str(batch) + '/' + f + '.csv', 'w')
        # fp_replicates.write('Drug, Concentration, Replicate_1, Replicate_2\n')

        # print f
        screen_results = get_feature_result(f, table)
        plates = screen_results.keys()

        # print screen_results

        # plate_means means of all DMSO wells on one plate
        drugs = {}
        for plate in plates:
            # plate well contains all the dmso (POC) values

            # go threw all the wells and look for DMSO treated wlls
            for well in screen_results[plate]:
                if screen_results[plate][well]['Worked'] == 'FALSE' or screen_results[plate][well][
                    'N_Value'] == -100000:
                    continue

                # if screen_results[plate][well]['Drug_2'] == 'DMSO':
                if screen_results[plate][well]['Drug_2'] == 'DMSO':
                    if drugs.has_key(screen_results[plate][well]['Drug_1']):
                        drugs[screen_results[plate][well]['Drug_1']].append(screen_results[plate][well]['N_Value'])
                    else:
                        drugs[screen_results[plate][well]['Drug_1']] = [screen_results[plate][well]['N_Value']]

        X = []
        Y = []
        for d, v in drugs.iteritems():

            drugs[d] = reject_outliers_2(drugs[d])

            if len(drugs[d]) <= 1:
                continue

            for i in range(0, len(drugs[d])):
                for i2 in range(0, len(drugs[d])):
                    if i > i2:
                        X.append(drugs[d][i])
                        Y.append(drugs[d][i2])

        # calculate pearson
        cor = np.corrcoef(X, Y)[0][1]
        # print f + ': %f' % cor
        fp.write(f+','+str(cor)+'\n')

        #In case a feature with only same values 'nan' is result/WIll be filtured through X_min = X_Max filter
        if str(cor) != 'nan':

            correlations.append(cor)

            # make plot
            if plot:
                plt.xlabel('Replicate 1')
                plt.ylabel('Replicate 2')
                plt.scatter(X, Y, s=2)
                plt.title(f)
                plt.legend(['Correlation = %.2f' % (cor)])
                #ensure_dir('../results/Correlation_Singles_Normalized/Correlation_' + f + '.pdf')
                # plt.savefig('../results/Correlation_Singles_Normalized/Correlation_' + f + '.pdf', format='pdf', dpi=800)
                plt.show()
                plt.close()



    plt.hist(correlations, bins='auto', color='grey')
    plt.axvline(0.3, ls='--', c='red')
    plt.legend(['Rejected Features: %d' % len([x for x in correlations if x < 0.3])])
    plt.savefig('../results/' + table + '/Correlation_Between_Replicates/Correlations.pdf')
    plt.close()

calculate_Correlation_BetweenPlates(table, False)