'''
This code includes all the information about the HCS screen and creates an input file that can be used for submission on the cluster


'''

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

def letter_to_Number(letter):
    #column = 1 to 24
    #row = A to P (01 to 16)

    dic_value = {'A':'01','B':'02','C':'03','D':'04','E':'05','F':'06','G':'07','H':'08','I':'09','J':'10','K':'11','L':'12','M':'13','N':'14','O':'15','P':'16' }
    return dic_value[letter]



if __name__ == "__main__":
    batch = 2


    plates_name = {}

    #load the input file
    fp_plate = open('../data/ImageInputFile/plateIDs_batch'+str(batch),'r')


    #get plate information
    for line in fp_plate:
        plates_name[line[0:7]] = line.strip()

    #check for correct batch
    if batch == 1:
        basis_path = '/scratch/lab_menche/CLOUD_imaging_001_to_064/'
    elif batch == 2:
        basis_path = '/scratch/lab_menche/CLOUD_imaging_065_to_124/'

    if batch == 1:
        fp = open('../data/ImageInputFile/20140206_transferlist_001_064_annotated.txt','r')
    elif batch == 2:
        fp = open('../data/ImageInputFile/20140225_transferlist_065_124_annotated.txt', 'r')

    #create directory if doesn't exist
    ensure_dir('../results/InputFile/InputFile_'+str(batch)+'.csv')
    fp_out = open('../results/InputFile/InputFile_'+str(batch)+'.csv','w')

    fp_out.write('URL_BetaTubulin,URL_DAPI,URL_MitoTracker,Metadata_Row,Metadata_Column,Metadata_Field,Metadata_Well,Metadata_Plate,Metadata_Transfer_A,Metadata_ID_A,Metadata_Name_A,Metadata_Transfer_B,Metadata_ID_B,Metadata_Name_B\n')


    #create output file
    group_index = 1
    fp.next()
    for line in fp:
        tmp =  line.split('\t')
        destWell = tmp[0]

        well_nums = letter_to_Number(destWell[0])+destWell[1:]

        imagestring = 'r'+well_nums[0:2]+'c'+well_nums[2:4]


        destPlate = tmp[3]


        if tmp[4] == 'Cpd':
            transferA = tmp[6]
            ID_A = tmp[8]
            Name_A = tmp[9]

            if tmp[13] == 'Cpd':
                transferB = tmp[15]
                ID_B = tmp[16]
                Name_B = tmp[17]
            else:
                transferB = tmp[15]
                ID_B = tmp[13]
                Name_B = tmp[13]
        else:
            transferA = tmp[6]
            ID_A = tmp[4]
            Name_A = tmp[4]

            ID_B = 'None'
            Name_B = 'None'
            transferB = 'None'

        for f in range(1,5):
            fp_out.write(
                         'file:'+basis_path+plates_name[destPlate]+'/Images/'+imagestring+'f0'+str(f)+'p01rc2-ch1sk1fk1fl1.tiff,'+
                         'file:' + basis_path + plates_name[destPlate] + '/Images/' + imagestring + 'f0' + str(f) + 'p01rc1-ch1sk1fk1fl1.tiff,' +
                         'file:' + basis_path + plates_name[destPlate] + '/Images/' + imagestring + 'f0' + str(f) + 'p01rc3-ch1sk1fk1fl1.tiff,' +
                         well_nums[0:2]+','+well_nums[2:4]+',' + '0'+str(f)+','+destWell+','+destPlate+','+transferA+','+ID_A+','+Name_A+','+transferB
                         + ',' + ID_B+','+Name_B+'\n')

            group_index+=1
    fp.close()
    fp_out.close()
