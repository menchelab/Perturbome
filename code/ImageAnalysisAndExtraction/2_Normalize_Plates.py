import MySQLdb
import pandas
import os

db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")

import numpy as np
from matplotlib import pylab as plt


class MedianPolish:
    """Fits an additive model using Tukey's median polish algorithm"""

    def __init__(self, array):
        """Get numeric data from numpy ndarray to self.tbl, keep the original copy in tbl_org"""
        if isinstance(array, np.ndarray):
            self.tbl_org = array
            self.tbl = self.tbl_org.copy()
        else:
            raise TypeError('Expected the argument to be a numpy.ndarray.')

    @staticmethod
    def csv_to_ndarray(fname):
        """ Utility method for loading ndarray from .csv file"""
        try:
            return np.genfromtxt(fname, delimiter=",")
        except Exception, e:
            print "Error loading file %s:" % fname
            raise

    def median_polish(self, max_iterations=10, method='median'):
        """
            Implements Tukey's median polish alghoritm for additive models
            method - default is median, alternative is mean. That would give us result equal ANOVA.
        """

        grand_effect = 0
        median_row_effects = 0
        median_col_effects = 0
        row_effects = np.zeros(shape=self.tbl.shape[0])
        col_effects = np.zeros(shape=self.tbl.shape[1])

        for i in range(max_iterations):
            if method == 'median':
                row_medians = np.median(self.tbl, 1)
                row_effects += row_medians
                median_row_effects = np.median(row_effects)

            elif method == 'average':
                row_medians = np.average(self.tbl, 1)
                row_effects += row_medians
                median_row_effects = np.average(row_effects)
            grand_effect += median_row_effects
            row_effects -= median_row_effects
            self.tbl -= row_medians[:, np.newaxis]

            if method == 'median':
                col_medians = np.median(self.tbl, 0)
                col_effects += col_medians
                median_col_effects = np.median(col_effects)
            elif method == 'average':
                col_medians = np.average(self.tbl, 0)
                col_effects += col_medians
                median_col_effects = np.average(col_effects)

            self.tbl -= col_medians

            grand_effect += median_col_effects

        return grand_effect, col_effects, row_effects, self.tbl, self.tbl_org


def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def loadBadImages(db_table):

    badImages = []
    fp = open('../results/'+db_table+'/BadImages/BadImages.csv', 'r')
    for line in fp:
        badImages.append(line.strip().split(',')[0])
    fp.close()

    return badImages


def Create_Table_Normlized_Well_Mean_FromImage_Table(db_table='Default'):
    '''
    This function normilzes wells according to the median subtraction methodology. This means that all wells get
    their corresponding meadian of it's row and column subtracted. This leads to normilzed values around 0 (for the median)
     > 0 if their value is bigger than the respective median or <0 if smaller. For values where the overall row/column
     shows higher values due to some plate effect, this subtraction is bigger leading to a removal of such plate effects.

     This function uses already pre calculated Mean values from cell_profiler (inside the Menche_DB)
    :return:
    '''

    show_PlateCorrection = True
    use_median_polish = True


    #Get BadImages
    badImages = loadBadImages(db_table)


    # Define table name
    if db_table == 'Default':
        print 'Set correct table:'
        exit()

    experiment_name = db_table.split('_')[0]

    # Get features
    string = "select COLUMN_NAME from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME='" + db_table + "'"
    all_features = list(pandas.read_sql(string, con=db)['COLUMN_NAME'])

    # Remove MetadataFeatures such as height of image
    features = []
    for f in all_features:
        # if 'Granularity' in f or 'Mean' in f or 'Median' in f:
        if ('Mean_' in f or 'Median_' in f) and 'Location' not in f and 'Center' not in f:
            features.append(f)
    features.sort()


    # Get all plates for this screen
    string = 'select Image_Metadata_Plate  from ' + db_table + ' group by Image_Metadata_Plate;'
    data_plates = pandas.read_sql(string, con=db)
    plates = data_plates['Image_Metadata_Plate']

    print 'Number of Features to calculate: %d' % len(features)
    print 'Number of Plates to normalize: %d' % len(plates)
    print 'Approximate time to complete: %.2f hours' % (0.5 * len(plates) * len(features) / 60 / 60)
    print '-------'

    # Calculate the individual normalized results for the features
    for f in features[0:1]:
        print f


        if os.path.isfile('../results/' + db_table + '/NormalizedFeatures/' + f + '.csv') == True:
            continue


        # Create Output File
        ensure_dir('../results/' + db_table + '/NormalizedFeatures/' + f + '.csv')
        fp_out = open('../results/' + db_table + '/NormalizedFeatures/'+ f + '.csv', 'w')
        fp_out.write('Image_Metadata_Plate,Image_Metadata_Well,Image_Metadata_ID_A,Image_Metadata_ID_B,Worked,' + f + '_Norm,' + f + '_MedianNorm,NotNormed,DMSO_Median,DMSO_MAD\n')


        #Get the AVG of one Well (ideally all 4 images of one well)

        f = 'AVG(' + f + ')'

        string = 'select COUNT(*),Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Transfer_A,Image_Metadata_Transfer_B,Image_Metadata_Plate,Image_Metadata_Well,' + f + ' from ' + db_table + ' where ImageNumber not in (' + ','.join(badImages) + ') group by Image_Metadata_Plate,Image_Metadata_Well;'
        data_plates = pandas.read_sql(string, con=db)


        # Go threw all plates
        for plate in plates:

            plate = 1315086
            data_well = data_plates.loc[data_plates['Image_Metadata_Plate'] == plate]

            # Get data for this plate
            #string = 'select Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Transfer_A,Image_Metadata_Transfer_B,Image_Metadata_Well,' + f + ' from ' + db_table + ' where Image_Metadata_Plate = ' + str(
            #    plate) + ' group by Image_Metadata_Well;'

            #string = 'select COUNT(*),Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Transfer_A,Image_Metadata_Transfer_B,Image_Metadata_Well,' + f + ' from ' + db_table + ' where Image_Metadata_Plate = ' + str(
            #    plate) + ' and ImageNumber not in (' + ','.join(badImages) + ') group by Image_Metadata_Well;'
            #data_well = pandas.read_sql(string, con=db)

            # Get all possible wells for this plate (sometimes not all 384 wells used)
            wells = list(set(data_well['Image_Metadata_Well']))
            wells.sort()

            print wells
            exit()

            # Create pandas dataframe, used for later normilization
            columns = []
            for i in range(1, 25):
                columns.append('%02.d' % i)
            df = pandas.DataFrame(columns=columns,
                                  index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                         'P'])

            dmso_wells = []
            # Fill pandas dataframe
            for well in wells:
                df[well[1:3]][well[0]] = data_well.loc[data_well['Image_Metadata_Well'] == well][f].values[0]
                # add dmso data to dmso list
                if data_well.loc[data_well['Image_Metadata_Well'] == well]['Image_Metadata_ID_A'].values[0] == 'DMSO':
                    dmso_wells.append(data_well.loc[data_well['Image_Metadata_Well'] == well][f].values[0])
            dmso_wells = [x for x in dmso_wells if str(x) !='nan']


            if use_median_polish:
                df = df.fillna(np.mean(dmso_wells))

                arr = np.array(df)

                # print arr
                mp = MedianPolish(arr)
                ge, ce, re, resid, tbl_org = mp.median_polish(20)


                if show_PlateCorrection:
                    print "median polish:"
                    print "grand effect = ", ge
                    print "column effects = ", ce
                    print "row effects = ", re
                    print "-----Table of Residuals-------"
                    #print resid

                    plt.title('Show Plate Correction')
                    plt.subplot(121)
                    plt.title('Corrected')
                    plt.imshow(resid, cmap='hot', interpolation='nearest')
                    print "-----Original Table-------"
                    #print tbl_org
                    plt.subplot(122)
                    plt.title('Original')
                    plt.imshow(tbl_org, cmap='hot', interpolation='nearest')
                    plt.show()
                    plt.close()


                df = pandas.DataFrame(data=resid, columns=['%02.d' % i for i in range(1, 25)],
                                      index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                             'P'])
                df_org = pandas.DataFrame(data=tbl_org, columns=['%02.d' % i for i in range(1, 25)],
                                          index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                                                 'O', 'P'])

            # perform normilzation for each well
            for well in wells:

                print well

                dmso_wells = [k for k in dmso_wells if str(k) != 'nan']
                DMSO_MAD = np.mean(np.abs(np.tile(np.mean(dmso_wells), (1, len(dmso_wells))) - dmso_wells))
                DMSO_median = np.median(dmso_wells)
                # DMSO_iqr = np.subtract(*np.percentile(dmso_wells, [75, 25]))
                # DMSO_Variation = DMSO_iqr / DMSO_MAD

                Median_normalized_value = (df_org[well[1:3]][well[0]] - DMSO_median) / (DMSO_MAD * 1.4826)

                # print well
                # print df
                # Get the value of the well
                current_value = df[well[1:3]][well[0]]

                if use_median_polish:
                    normalized_value = current_value
                    current_value = df_org[well[1:3]][well[0]]
                    values = dmso_wells

                else:
                    # Get the value of the well
                    #current_value = df[well[1:3]][well[0]]

                    # get all values of column and rows, remove 'nan'
                    values = list(df.loc[well[0]])
                    values.extend(df[well[1:3]])
                    values = [k for k in values if str(k) != 'nan']

                    # Subtract Median from Well value
                    normalized_value = (current_value - np.median(values))

                # Extract further information
                Trans_A = list(data_well.loc[data_well['Image_Metadata_Well'] == well]['Image_Metadata_Transfer_A'])[0]
                Trans_B = list(data_well.loc[data_well['Image_Metadata_Well'] == well]['Image_Metadata_Transfer_B'])[0]

                ID_A = list(data_well.loc[data_well['Image_Metadata_Well'] == well]['Image_Metadata_ID_A'])[0]
                ID_B = list(data_well.loc[data_well['Image_Metadata_Well'] == well]['Image_Metadata_ID_B'])[0]

                count = list(data_well.loc[data_well['Image_Metadata_Well'] == well]['COUNT(*)'])[0]

                if count <= 2:
                    worked = 'FALSE'
                elif ID_A == 'DMSO' or ID_A == 'PosCon':
                    if Trans_A == 'YES':
                        worked = 'TRUE'
                    else:
                        worked = 'FALSE'
                else:
                    if Trans_A == 'YES' and Trans_B == 'YES':
                        worked = 'TRUE'
                    else:
                        worked = 'FALSE'

                if str(normalized_value) == 'nan':
                    normalized_value = -100000
                    worked = 'FALSE'

                # Write Output
                fp_out.write(str(plate) + ',' + well + ',' + ID_A + ',' + ID_B + ',' + worked + ',' + str(
                    normalized_value) + ',' + str(Median_normalized_value) + ',' + str(current_value) + ',' + str(
                    np.median(values)) + ',' + str(DMSO_MAD) + '\n')


if __name__ == "__main__":
    Create_Table_Normlized_Well_Mean_FromImage_Table('DPN1018Batch2Per_Image')



    #Median_Cells_Texture_SumVariance_DAPI_3_03.csv
    #Median_Nuclei_AreaShape_Zernike_4_4.csv