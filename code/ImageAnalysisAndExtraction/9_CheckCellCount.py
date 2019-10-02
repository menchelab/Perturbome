import pandas
import MySQLdb
import numpy as np
import pickle
import os

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def create_Single_CellCounts(db_table):
    db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")

    string = "select Image_Metadata_ID_A from DDI0517Per_Image group by Image_Metadata_ID_A;"

    data = pandas.read_sql(string, con=db)['Image_Metadata_ID_A']

    #with open('../results/FeatureVectors/SingleVectors_' + str(min(plates)) + '_to_' + str(
    #        max(plates)) + '_NoCutoff_' + str(cast_int) + '.pickle', 'rb') as handle:
    #    single_Vectors = pickle.load(handle)

    singles = list(data)
    singles.sort()

    if 'PosCon' in singles:
        singles.remove('PosCon')

    if 'DMSO' in singles:
        singles.remove('DMSO')

    # Define Database to check for missing Images

    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Plate from " + db_table + " where Image_Metadata_ID_B like 'DMSO'  and Image_Metadata_Transfer_A  like 'YES' and Image_Metadata_Transfer_B like 'YES' group by Image_Metadata_ID_A,Image_Metadata_Plate,Image_Metadata_Well;"

    data = pandas.read_sql(string,con=db)

    ensure_dir('../results/'+table+'/CellCount/SinglesCellCount.csv')
    fp_out = open('../results/'+table+'/CellCount/SinglesCellCount.csv','w')
    fp_out.write('Drug,AVG_CellCount\n')
    for drug in singles:

        if len(data.loc[data['Image_Metadata_ID_A'] == drug]['SUM(Image_Count_Cells)'].values) > 0:
            cellcount =  np.mean(data.loc[data['Image_Metadata_ID_A'] == drug]['SUM(Image_Count_Cells)'].values)
            cellcount = int(cellcount)
        else:
            cellcount = 'nan'

        fp_out.write(drug+','+str(cellcount) +'\n')
    fp_out.close()

def create_Single_CellCounts_individualReplicates(db_table):
    db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")

    string = "select Image_Metadata_ID_A from DDI0517Per_Image group by Image_Metadata_ID_A;"

    data = pandas.read_sql(string, con=db)['Image_Metadata_ID_A']

    #with open('../results/FeatureVectors/SingleVectors_' + str(min(plates)) + '_to_' + str(
    #        max(plates)) + '_NoCutoff_' + str(cast_int) + '.pickle', 'rb') as handle:
    #    single_Vectors = pickle.load(handle)

    singles = list(data)
    singles.sort()

    if 'PosCon' in singles:
        singles.remove('PosCon')

    if 'DMSO' in singles:
        singles.remove('DMSO')


    plates = range(1315001, 1315124, 10)

    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Plate from " + db_table + " where Image_Metadata_ID_B like 'DMSO'  and Image_Metadata_Transfer_A  like 'YES' and Image_Metadata_Transfer_B like 'YES' group by Image_Metadata_ID_A,Image_Metadata_Plate,Image_Metadata_Well;"

    data = pandas.read_sql(string,con=db)

    ensure_dir('../results/' + table + '/CellCount/SinglesCellCount_AllReplicates.csv')
    fp_out = open('../results/' + table + '/CellCount/SinglesCellCount_AllReplicates.csv','w')
    #fp_out.write('Drug,CellCounts\n')
    fp_out.write('Drug,' +','.join([str(x) for x in plates]) +'\n')
    for drug in singles:

        cellcounts = []
        for p in plates:
            value =  data.loc[(data['Image_Metadata_ID_A'] == drug) & (data['Image_Metadata_Plate'] == p)]['SUM(Image_Count_Cells)'].values
            if len(value) != 0:
                cellcounts.append(str(value[0]))
            else:
                cellcounts.append('None')

        fp_out.write(drug + ',' + ','.join([str(x) for x in cellcounts]) + '\n')


    fp_out.close()

def create_Combination_CellCounts(db_table):
    # Define Database to check for missing Images
    db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")

    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Plate,Image_Metadata_Transfer_A,Image_Metadata_Transfer_B from " + db_table + " where Image_Metadata_ID_B like 'CLOUD%' group by Image_Metadata_ID_A,Image_Metadata_ID_B,Image_Metadata_Well,Image_Metadata_Plate;"
    data = pandas.read_sql(string,con=db)

    ensure_dir('../results/' + table + '/CellCount/CombisCellCount.csv')
    fp_out = open('../results/' + table + '/CellCount/CombisCellCount.csv','w')
    fp_out.write('Plate,Well,Drug1,Drug2,AVG_CellCount\n')
    for c in data.iterrows():

        if c[1][5] == 'YES' and c[1][6] == 'YES':
            value = str(c[1][0])
        else:
            value = 'nan'

        fp_out.write(str(c[1][4])+','+c[1][1]+',' +c[1][2]+','+c[1][3]+','+value +'\n')
    fp_out.close()

def getDMSO_Untreated_CellCount(db_table):
    # Define Database to check for missing Images
    db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")


    string = "select SUM(Image_Count_Cells), Image_Metadata_Well, Image_Metadata_Plate from " + db_table + "  where Image_Metadata_ID_A like 'DMSO' and Image_Metadata_ID_B like 'None' and Image_Metadata_Transfer_A  like 'YES' group by Image_Metadata_Well,Image_Metadata_Plate;"
    data = pandas.read_sql(string,con=db)



    mean = np.mean(data['SUM(Image_Count_Cells)'])
    std =  np.std(data['SUM(Image_Count_Cells)'])
    max_val = np.percentile(data['SUM(Image_Count_Cells)'],98)

    ensure_dir('../results/' + table + '/CellCount/DMSO_Overview.csv')
    fp_out = open('../results/' + table + '/CellCount/DMSO_Overview.csv', 'w')
    fp_out.write('Mean,Std,Max\n%f,%f,%f' %(mean,std,max_val))
    fp_out.close()

    fp_out = open('../results/' + table + '/CellCount/DMSO_Replicates.csv', 'w')
    fp_out.write('Plate,Well,CellCount\n')
    for row in data.iterrows():
        fp_out.write(str(row[1][2])+','+row[1][1]+','+str(row[1][0])+'\n')
    fp_out.close()

def get_CellCount_perWell(db_table):
    # Define Database to check for missing Images
    db = MySQLdb.connect("menchelabdb.int.cemm.at", "root", "cqsr4h", "ImageAnalysisDDI")


    string = "select SUM(Image_Count_Cells),Image_Metadata_ID_A,Image_Metadata_ID_B, Image_Metadata_Well, Image_Metadata_Plate,Image_Metadata_Transfer_A,Image_Metadata_Transfer_B from " + db_table + " group by Image_Metadata_Well,Image_Metadata_Plate;"
    data = pandas.read_sql(string,con=db)

    data.sort_values(by=['Image_Metadata_Plate','Image_Metadata_Well'])

    ensure_dir('../results/' + db_table + '/CellCount/Individual_Well_Results.csv')
    fp_out = open('../results/' + db_table + '/CellCount/Individual_Well_Results.csv', 'w')

    fp_out.write('ID_A,ID_B,Plate,Well,CellCount,TransferOK\n')
    for row in data.iterrows():

        ID_A = row[1][1]
        ID_B = row[1][2]

        Trans_A = row[1][5]
        Trans_B = row[1][5]


        if ID_A == 'DMSO' or ID_A == 'PosCon':
            if Trans_A == 'YES':
                worked = 'TRUE'
            else:
                worked = 'FALSE'
        else:
            if Trans_A == 'YES' and Trans_B == 'YES':
                worked = 'TRUE'
            else:
                worked = 'FALSE'


        fp_out.write(ID_A+','+ID_B+','+str(row[1][4])+','+row[1][3]+','+str(row[1][0])+','+worked+'\n')
    fp_out.close()



def Combine_Well_Results(table1,table2):

    if os.path.isfile('../results/' + table1 + '/CellCount/Individual_Well_Results.csv') == False or os.path.isfile('../results/' + table2 + '/CellCount/Individual_Well_Results.csv') == False:
        print 'First create single files'
        exit()
    else:

        fp_out = open('../results/All_CellCounts_Combined.csv', 'w')

        fp = open('../results/' + table1 + '/CellCount/Individual_Well_Results.csv')
        for line in fp:
            fp_out.write(line)
        fp.close()

        fp = open('../results/' + table2 + '/CellCount/Individual_Well_Results.csv')
        fp.next()
        for line in fp:
            fp_out.write(line)
        fp.close()

        fp_out.close()

def PlotResult_file(all=False):


    from matplotlib import pylab as plt

    drug_values = {}
    dmso_values = []
    fp = open('../results/' + table1 + '/CellCount/Individual_Well_Results.csv')
    fp.next()
    for line in fp:


        tmp = line.strip().split(',')
        if tmp[5] == 'TRUE':

            if all or tmp[1] == 'DMSO':

                if drug_values.has_key(tmp[0]):
                    drug_values[tmp[0]].append(float(tmp[4]))
                else:
                    drug_values[tmp[0]] = [float(tmp[4])]

            if tmp[0] == 'DMSO':
                dmso_values.append(float(tmp[4]))

    max_val  = np.mean(dmso_values) + 0.5 * np.std(dmso_values)
    #max_val =  np.mean([np.mean(x) for x in drug_values.values()]) + 1.2 * np.std([np.mean(x) for x in drug_values.values()])


    effect = 0

    normalized = []
    for drug in drug_values:
        scaled = (np.mean(drug_values[drug]) - 0) / max_val
        if scaled <= 1:
            normalized.append(scaled)
        else:
            normalized.append(1)

        if scaled < 0.5:
            effect +=1

    print 'Number of drugs with more than 50%% cytotoxicity: %d' %effect
    print  'Number of drugs with les  than 50%% cytotoxicity: %d' %(len(drug_values) - effect)

    plt.hist(normalized,bins='auto', color = '#40B9D4')
    #plt.show()
    plt.savefig('../results/CellCountHistogram.pdf')
    plt.close()

# - - - - - - - - - - - - - - - - - - - -
# Define Experiment
table1 = 'DPN1018Batch1Per_Image'
table2 = 'DPN1018Batch2Per_Image'
# - - - - - - - - - - - - - - - - - - - -



#create_Single_CellCounts(table)
#create_Single_CellCounts_individualReplicates(table)
#getDMSO_Untreated_CellCount(table)
#create_Combination_CellCounts(table)

#get_CellCount_perWell(table1)

#Combine_Well_Results(table1,table2)


PlotResult_file()