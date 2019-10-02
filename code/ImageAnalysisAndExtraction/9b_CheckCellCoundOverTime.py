import pandas
import MySQLdb
from matplotlib import pylab as plt
import numpy as np
import os

db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)



def Check_Drug_Decay(db_table):
    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Plate from " + db_table + " where Image_Metadata_ID_B like 'DMSO'  group by Image_Metadata_ID_A,Image_Metadata_Plate,Image_Metadata_Well;"
    data = pandas.read_sql(string, con=db)

    all_clouds = list(set(data['Image_Metadata_ID_A']))
    all_clouds.sort()


    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Plate from " + db_table + " where Image_Metadata_ID_B like 'DMSO'  and Image_Metadata_Transfer_A  like 'YES' and Image_Metadata_Transfer_B like 'YES' group by Image_Metadata_ID_A,Image_Metadata_Plate,Image_Metadata_Well;"
    data = pandas.read_sql(string, con=db)


    # plates = timepoints
    plates = list(set(data['Image_Metadata_Plate']))
    plates.sort()

    ensure_dir('../results/' + db_table + '/DrugDecay_OverTime/Overview.csv')
    fp_out = open('../results/' + db_table + '/DrugDecay_OverTime/Overview.csv', 'w')
    fp_out.write('Drug,Steep,MaxDifference\n')

    max_val = np.mean(data['SUM(Image_Count_Cells)']) + 1.95 * np.std(data['SUM(Image_Count_Cells)'])

    for c in all_clouds:

        print c
        drug_values = data.loc[data['Image_Metadata_ID_A'] == c]

        x = []
        y = []
        for timepoint, p in enumerate(plates):

            val = drug_values.loc[drug_values['Image_Metadata_Plate'] == p]['SUM(Image_Count_Cells)'].values
            if len(val) > 0:
                scaled = (val[0] - 0) / max_val
                x.append(scaled)
                y.append(timepoint)


        if len(x) > 1:
            biggest_difference = max(x) - min(x)
        else:
            biggest_difference = 'nan'

        if len(x) < 3:
            fp_out.write(c + ',nan,'+str(biggest_difference)+'\n')
        else:

            fit = np.polyfit(y, x, 1)



            y1 = 0 * (fit[0]) + fit[1]
            y2 = (len(y) - 1) * (fit[0]) + fit[1]

            fp_out.write(c + ',' + str(fit[0])+','+str(biggest_difference) + '\n')

        if plot and abs(fit[0]) < 0.05:
            plt.plot([0, len(y) - 1], [y1, y2], c='#40B9D4')
            plt.legend(['y = %.2fx+ %.2f' % (fit[0], fit[1])])
            plt.scatter(y, x, c='grey')

            plt.title(c)
            plt.xlabel('Timepoints')
            plt.ylabel('Cell Count')
            plt.ylim(0, 1)

            # ensure_dir()

            plt.savefig('../results/' + db_table + '/DrugDecay_OverTime/' + c + '.pdf', dpi=600)
            # plt.show()
            plt.close()

    fp_out.close()


def Combine_Well_Results(table1,table2):

    if os.path.isfile('../results/' + table1 + '/DrugDecay_OverTime/Overview.csv') == False or os.path.isfile('../results/' + table2 + '/DrugDecay_OverTime/Overview.csv') == False:
        print 'First create single files'
        exit()
    else:



        drug_decays = {}
        drug_maxDifference = {}

        drugs = set()

        drug_decays[table1] = {}
        drug_maxDifference[table1] = {}

        fp = open('../results/' + table1 + '/DrugDecay_OverTime/Overview.csv')
        fp.next()
        for line in fp:
            tmp = line.strip().split(',')

            print tmp
            drugs.add(tmp[0])
            drug_decays[table1][tmp[0]] = tmp[1]
            drug_maxDifference[table1][tmp[0]] = tmp[2]

        fp.close()

        drug_decays[table2] = {}
        drug_maxDifference[table2] = {}

        fp = open('../results/' + table2 + '/DrugDecay_OverTime/Overview.csv')
        fp.next()
        for line in fp:
            tmp = line.strip().split(',')
            drugs.add(tmp[0])
            drug_decays[table2][tmp[0]] = tmp[1]
            drug_maxDifference[table2][tmp[0]] = tmp[2]

        fp.close()

        fp_out = open('../results/DrugDecay_Combined.csv', 'w')
        fp_out.write('Drug,Batch1_Decay,Batch1_Diff,Batch2_Decay,Batch2_Diff\n')
        for drug in drugs:

            decay_1 = 'nan'
            if drug_decays[table1].has_key(drug):
                decay_1 = drug_decays[table1][drug]

            decay_2 = 'nan'
            if drug_decays[table2].has_key(drug):
                decay_2 = drug_decays[table2][drug]

            diff_1 = 'nan'
            if drug_maxDifference[table1].has_key(drug):
                diff_1 = drug_maxDifference[table1][drug]

            diff_2 = 'nan'
            if drug_maxDifference[table2].has_key(drug):
                diff_2 = drug_maxDifference[table2][drug]

            fp_out.write(drug+','+decay_1+','+diff_1+','+decay_2+','+diff_2+'\n')
        fp_out.close()


# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
db_table1 = 'DPN1018Batch1Per_Image'
db_table2 = 'DPN1018Batch2Per_Image'
plot = True
# - - - - - - - - - - - - - - - - - - - -

Check_Drug_Decay(db_table1)

#Combine_Well_Results(db_table1,db_table2)