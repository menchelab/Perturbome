'''
THIS SCRIPT FINDS MISSING IMAGES IN THE DATABASES
THIS COULD OCCUR BECAUSE THE CLUSTER ANALYSIS MIGHT ENCOUNTERED SOME BUG OR RAN OUT OF TIME FOR ANALYSIS

SIMPLE ALGORITHM THAT SEARCHES FOR HOLES IN IMAGE NUMBERS
E.g. 1,2,3,4,8,9,10,11... --> Missing Images between 5 and 8 (BOTH NUMBER INCLUDED AS MISSING!)
'''


import pandas
import MySQLdb
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

#Define Database to check for missing Images
db = MySQLdb.connect("menchelabdb.int.cemm.at","root","cqsr4h","ImageAnalysisDDI" )

if __name__ == "__main__":
    #Define Table for checking for ImageNumber
    string = "select ImageNumber from DPN1018Batch2Per_Image;"

    #Load Image Numbers
    data = pandas.read_sql(string, con=db)

    #Create Directory if doesn't exist
    ensure_dir('../results/ImageInputFile/missing_images.csv')

    #Open output file
    fp = open('../results/ImageInputFile/missing_images.csv','w')

    #Iterate threw Image Numbers and check for holes (= if current != line)
    #if current is not line then the numbers in between are missing
    #set current to line
    #continue
    current = 1
    for line in data.ImageNumber:
        if line != current:
            #print 'Found Missing:'
            #print 'From: %d to %d' %(current,line-1)
            if (line-1) - current < 30:
                print 'sbatch UseInputFile.sh %d %d' %(current,line-1)
            else:
                split =  int(((line-1) - current)/25)
                for i in range(1,split):
                    print  'sbatch UseInputFile.sh %d %d' %(current,current+25)
                    current = current + 25 +1

                remaining = (line - 1) - current
                if remaining!=0:
                    print  'sbatch UseInputFile.sh %d %d' % (current, line - 1)



            fp.write(str(current)+','+str(line-1)+'\n')
            current = line

        current +=1

    fp.close()


    # SBATCH --mail-type=end
    # SBATCH --mail-user=mcaldera@cemm.oeaw.ac.at