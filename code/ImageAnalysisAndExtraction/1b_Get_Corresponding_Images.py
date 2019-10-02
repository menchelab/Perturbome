import pandas
import MySQLdb
from PIL import Image
from matplotlib import pylab as plt
import numpy as np
from scipy import misc
from skimage import exposure
from skimage.exposure import rescale_intensity
from matplotlib.backends.backend_pdf import PdfPages
import os

#Define Database
db = MySQLdb.connect("menchelabdb.int.cemm.at","root","cqsr4h","ImageAnalysisDDI" )

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def find_ImageByNumber(imageNumber,table):
    data = pandas.read_sql("select Image_FileName_DAPI,Image_FileName_Mitotracker,Image_FileName_BetaTubulin,Image_PathName_DAPI from "+table+" where ImageNumber  = '" + imageNumber + "';",con=db).values[0]


    path = data[3].split('lab_menche')[1]



    im = Image.open(
        "/Volumes/scratch/lab_menche" + path + "/" + data[0])
    im2 = Image.open(
        "/Volumes/scratch/lab_menche" + path + "/" + data[2])
    im3 = Image.open(
        "/Volumes/scratch/lab_menche" + path + "/" + data[1])

    k = np.array(im).astype(np.int16)
    k2 = np.array(im2).astype(np.int16)
    k3 = np.array(im3).astype(np.int16)

    rgb = np.dstack((k3, k2, k))
    r_rgb = misc.imresize(rgb, 100)

    test = misc.bytescale(r_rgb)

    plt.imshow(test)
    #plt.show()
    ensure_dir('../results/'+table+'/BadImages/Images/Image_'+str(imageNumber)+'.png')
    plt.savefig('../results/'+table+'/BadImages/Images/Image_'+str(imageNumber)+'.png')
    plt.close()


if __name__ == "__main__":

    table = 'DPN1018Batch1Per_Image'


    fp = open('../results/'+table+'/BadImages/BadImages.csv', 'r')
    for line in fp:
        tmp = line.strip().split(',')
        find_ImageByNumber(tmp[0], table)
    fp.close()
