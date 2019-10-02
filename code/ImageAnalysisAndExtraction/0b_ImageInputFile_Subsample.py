'''
Keep only a certain subsample of the input file for creating the pipeline

Input Files NEED to be first created
'''
import random
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

batch = 1

#Load the interesting drugs
drugs_of_interest =[]
fp = open('../data/Drugs_Of_Interest_Batch'+str(batch))
for line in fp:
    drugs_of_interest.append(line.strip())
fp.close()

ensure_dir('../results/InputFile/InputFile_SubSampled_Batch' + str(batch) + '.csv')
fp_out = open('../results/InputFile/InputFile_SubSampled_Batch'+str(batch)+'.csv','w')


#These results in ca. 100 Entries
ratio_keep_single = 150 #150 is good
ratio_keep_combination = 1 #1 is god
ratio_keep_DMSO = 2 #2 is good

fp = open('../results/InputFile_'+str(batch)+'.csv')
header = fp.readline()
fp_out.write(header)


prev_drug = ''
for line in fp:
    tmp =  line.strip().split(',')


    rand_num = random.randint(0,1000)

    if tmp[9] == 'DMSO':
        if rand_num < ratio_keep_DMSO:
            fp_out.write(line)
    elif tmp[9] in drugs_of_interest and tmp[12] == 'DMSO' and tmp[8] != prev_drug:
        if rand_num < ratio_keep_single:
            fp_out.write(line)
            prev_drug = tmp[9]
    elif tmp[9] in drugs_of_interest or tmp[12] in drugs_of_interest:
        if rand_num < ratio_keep_combination:
            fp_out.write(line)



fp.close()
fp_out.close()

    #print tmp[8]
    #print tmp[11]
