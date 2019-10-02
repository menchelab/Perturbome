import os
import numpy as np

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

#Methods for getting feature by feature
def getFeatureList(db_table):
    mypath = '../results/' + db_table + '/NormalizedFeatures/'

    #mypath = '../data/Normalized_Wells/'
    onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]      # get all files of this folder (but no directories)
    onlyfiles.sort()
    features = []
    for file in onlyfiles:
        features.append(file.strip().split('.')[0])

    if '' in features:
        features.remove('')

    return  features

#Get the actual feature results from the Tukey plate normalized values
def get_feature_result(feature,db_table):


    #go threw the input file of the feature
    #path = '../data/Normalized_Wells/' + feature + '.csv'
    path = '../results/'+db_table+'/NormalizedFeatures/' + feature + '.csv'
    fp = open(path, 'r')
    fp.next()
    feature_results = {}
    #mean = {}
    for line in fp:
        tmp = line.strip().split(',')


        plate = int(tmp[0])
        well = tmp[1]
        drug1 = tmp[2]
        drug2 = tmp[3]
        worked = tmp[4]


        #if 'nan' for some features this might happen, then just set to mean of the plate
        if tmp[5] != 'nan':
            normed_value = tmp[5]
        else:
            #normed_value = np.mean(mean[plate])
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


    return  feature_results

# Extract the maximum and minimum value for the individual features
def get_Xmin_Xmax_perFeatures(db_table):
    '''
    for each feature find the borders to be set.
    Right now the top and lower border are 0.5 percentile as well as 99.5. The most extreme percent (in both direction)
    might be an outlier.

    :param batch:
    :return:
    '''

    # load all features
    features = getFeatureList(db_table)

    #contains all the feature max/min
    feature_MaxMin = {}

    #output file: create a file for each batch that contains all Max and Min values per feature
    ensure_dir('../results/'+db_table+'/POCNormalized/MaxMin_Values.csv')
    fp = open('../results/'+db_table+'/POCNormalized/MaxMin_Values.csv','w')
    fp.write('Feature,Max,Min\n')

    # go through all the features
    for f in features:
        feature_MaxMin[f] = {'Max':0.0,'Min':0.0}
        # load the data for the spcific feature
        screen_results = get_feature_result(f,db_table)

        plates = screen_results.keys()

        #get all values
        all_values = []
        for plate in plates:
            all_values.extend([screen_results[plate][k]['N_Value'] for k in screen_results[plate] if
                          str(screen_results[plate][k]['N_Value']) != '-100000' and screen_results[plate][k]['Worked'] != 'FALSE'])

        if len(all_values) == 0:
            print 'No Valid results for: %s' %f
            continue
        max_val = np.percentile(all_values, 99.5)   # discard 0.5% of the highest & lowest values to avoid normalization according to outliers
        min_val = np.percentile(all_values, 0.5)


        feature_MaxMin[f]['Max'] = max_val
        feature_MaxMin[f]['Min'] = min_val
        fp.write(f+','+str(max_val)+','+str(min_val)+'\n')

    return feature_MaxMin

def Normalize_To_POC(db_table='Default'):

    # Define table name
    if db_table == 'Default':
        print 'Set correct table:'
        exit()

    if os.path.isfile('../results/'+db_table+'/POCNormalized/MaxMin_Values.csv') == False:
        print 'Extracting Min/Max Values'
        feature_MaxMin = get_Xmin_Xmax_perFeatures(db_table)
    else:
        feature_MaxMin = {}
        fp = open('../results/'+db_table+'/POCNormalized/MaxMin_Values.csv', 'r')
        fp.next()
        for line in fp:
            tmp = line.strip().split(',')
            feature_MaxMin[tmp[0]] = {'Max':float(tmp[1]),'Min':float(tmp[2])}


    for feature, values in feature_MaxMin.iteritems():
        # pick min & max values from MinMax_Values file
        max = values['Max']
        min = values['Min']


        if min != max:

            fp_POC = open('../results/'+db_table+'/POCNormalized/' + feature + '.csv','w')
            fp_POC.write('Image_Metadata_Plate,Image_Metadata_Well,Image_Metadata_ID_A,Image_Metadata_ID_B,Worked,POC_Normalized\n')

            # open feature file
            fp_features = open('../results/'+db_table+'/NormalizedFeatures/' + feature + '.csv', 'r')
            fp_features.next()


            for line_f in fp_features:
                tmp_f = line_f.strip().split(',')
                first_part = tmp_f[0:5]     # new/normalized value should be written between first and second part

                value_to_normalize = float(tmp_f[5])

                if value_to_normalize > max:    # all values that are higher than the max value should be set to 1
                    POC_normalized = 1
                elif value_to_normalize < min:      # all values that are smaller than the min value should be set to 0
                    POC_normalized = 0
                else:                               # if values are between max and min, they should be normalized between 0 and 1
                    POC_normalized = (value_to_normalize - min) / (max - min)

                fp_POC.write(','.join(first_part)+','+str(POC_normalized)+'\n')


# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
table = 'DPN1018Batch1Per_Image'

# - - - - - - - - - - - - - - - - - - - -



if __name__ == "__main__":
    print 'Normalizing to POC for: %s' %table

    Normalize_To_POC(table)
    #exit()
