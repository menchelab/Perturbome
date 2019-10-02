import pandas
import MySQLdb
from PIL import Image
from matplotlib import pylab as plt
import numpy as np
from scipy import misc
from skimage import exposure,img_as_int
from skimage.transform import resize

from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
import os



#Define Database
db = MySQLdb.connect("menchelabdb.int.cemm.at","root","cqsr4h","ImageAnalysisDDI" )


def create_combined_Image(images,foci):
    for i in range(0,len(images)):
        v_min, v_max = np.percentile(images[i].flatten(), (0.01, 99.99))
        better_contrast = exposure.rescale_intensity(images[i],in_range=(v_min, v_max))
        images[i] = Image.fromarray(better_contrast,mode='RGB')


    new_im = Image.new('RGB', (2000, 1500))
    coordinates = {2: (0, 0), 3: (1000, 0), 1: (0, 750), 4: (1000, 750)}
    for image, c in zip(images, foci):
        image.thumbnail((1000, 1000))
        new_im.paste(image, coordinates[c])

    #plt.imshow(new_im)
    #plt.show()

    return new_im

def find_Single(drug1,batch):


    data = pandas.read_sql("select Image_FileName_DAPI,Image_FileName_Mitotracker,Image_FileName_BetaTubulin,Image_PathName_DAPI  from DPN1018"+batch+"Per_Image where (Image_Metadata_ID_A  = '" + drug1 + "' and Image_Metadata_ID_B = 'DMSO');",con=db)


    if len(data) == 0:
        print 'No Images found (Check spelling)'
        exit()

    images = []
    foci = []
    for image, image2, image3,path in zip(data['Image_FileName_DAPI'], data['Image_FileName_BetaTubulin'], data['Image_FileName_Mitotracker'],data['Image_PathName_DAPI']):

        path = path.split('lab_menche')[1]

        foci.append(int(image.split('f0')[1].split('p')[0]))


        k  = misc.imread("/Volumes/scratch/lab_menche"+path+"/" + image)
        k2 = misc.imread("/Volumes/scratch/lab_menche" + path + "/" + image2)
        k3 = misc.imread("/Volumes/scratch/lab_menche" + path + "/" + image3)



        rgb = np.dstack((k3, k2, k))
        #r_rgb = misc.imresize(rgb, 100)

        test = misc.bytescale(rgb)

        #plt.imshow(test)
        #plt.show()
        #plt.close()

        images.append(test)



    combined_images = []
    for i in range(0,len(images),4):
        c = create_combined_Image(images[i:i+4], foci[i:i+4])
        combined_images.append(c)




    return combined_images

def find_Combination(drug1,drug2,batch):


    data = pandas.read_sql("select  Image_FileName_DAPI,Image_FileName_Mitotracker,Image_FileName_BetaTubulin,Image_PathName_DAPI  from DPN1018"+batch+"Per_Image where (Image_Metadata_ID_A  = '"+drug1+"' and Image_Metadata_ID_B = '"+drug2+"') or (Image_Metadata_ID_A  = '"+drug2+"' and Image_Metadata_ID_B = '"+drug1+"');", con=db)


    if len(data) == 0:
        print 'No Images found (Check spelling)'
        exit()

    images = []
    foci = []
    for image, image2,image3,path in zip(data['Image_FileName_DAPI'], data['Image_FileName_BetaTubulin'], data['Image_FileName_Mitotracker'],data['Image_PathName_DAPI']):
        path = path.split('lab_menche')[1]

        foci.append(int(image.split('f0')[1].split('p')[0]))

        k = misc.imread("/Volumes/scratch/lab_menche" + path + "/" + image)
        k2 = misc.imread("/Volumes/scratch/lab_menche" + path + "/" + image2)
        k3 = misc.imread("/Volumes/scratch/lab_menche" + path + "/" + image3)

        #print k
        #exit()

        rgb = np.dstack((k3, k2, k))



        #r_rgb = misc.imresize(rgb, 100)
        #plt.imshow(rgb)
        #plt.show()

        test = misc.bytescale(rgb)

        #plt.imshow(test)
        #plt.show()
        #plt.close()

        images.append(test)




    return create_combined_Image(images,foci)



'''
GET ALL COMBINATIONS THAT SHOULD GET IMAGES RETRIEVED
'''

#combi[batch]
combinations = {}
fp = open('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/InteractionsDetailedOverview.csv','r')
fp.next()
for line in fp:
    tmp = line.strip().split(',')
    combinations[(tmp[0]),tmp[1]] = {'Batch':tmp[2],'Type':tmp[3],'Strength':tmp[4]}




alreadyExisting = [f for f in os.listdir('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/Interaction_Images/') if os.path.isfile(os.path.join('../results/Create_DPI_Network/Interaction_Images/', f))]


already_found = {}
for comb in list(combinations.keys()):

    try:
        d1 = comb[0]
        d2 = comb[1]

        batch = combinations[comb]['Batch']
        IntType =  combinations[comb]['Type']
        Strength = combinations[comb]['Strength']

        d1 = 'CLOUD057'

        if  d1 + ',' + d2 + '.pdf' in alreadyExisting:
            continue



        if d1 > d2:
            pca_plot_name = d1+'_'+d2+'.png'
        else:
            pca_plot_name = d2 + '_' + d1 + '.png'


        if os.path.isfile('../results/Calculate_Interactions/Combinations/'+pca_plot_name):
            PCA_Plot =   misc.imread('../results/Calculate_Interactions/Combinations/'+pca_plot_name, mode='RGB')
            PCA_Plot = Image.fromarray(PCA_Plot, mode='RGB')
            PCA_Plot = PCA_Plot.resize((1000,750))
    
    
            combination_image = find_Combination(d1, d2, batch)
            combination_images = Image.new('RGB', (1010 * 2 - 10, 750))
    
            combination_image.thumbnail((1000, 1000))
            combination_images.paste(combination_image, (0 * 1010, 0))
            PCA_Plot.thumbnail((1000, 1000))
            combination_images.paste(PCA_Plot, (1 * 1010, 0))
        else:

            combination_images = find_Combination(d1, d2, batch)


        if already_found.has_key(d1+batch) == True:
            d1_images = already_found[d1+batch]
        else:
            d1_images = find_Single(d1, batch)
            already_found[d1+batch] = d1_images

        if already_found.has_key(d2+batch) == True:
            d2_images = already_found[d2+batch]
        else:
            d2_images = find_Single(d2, batch)
            already_found[d2+batch] = d2_images

        combination_image = find_Combination(d1,d2,batch)

        s1_images = d1_images
        drug1_images = Image.new('RGB', (1010 * len(s1_images) - 10, 750))
        for i, image in enumerate(s1_images):
            image.thumbnail((1000, 1000))
            drug1_images.paste(image, (i * 1010, 0))

        s2_images = d2_images
        drug2_images = Image.new('RGB', (1010 * len(s2_images) - 10, 750))
        for i, image in enumerate(s2_images):
            image.thumbnail((1000, 1000))
            drug2_images.paste(image, (i * 1010, 0))

        plt.subplot(311)
        plt.imshow(combination_images)
        plt.ylabel('Combination (%s)' %batch)
        plt.xlabel(IntType +'(Maha_Dist: ' + Strength +')', fontsize=6)
        plt.subplot(312)
        plt.imshow(drug1_images)
        plt.ylabel(d1)
        plt.subplot(313)
        plt.imshow(drug2_images)
        plt.ylabel(d2)
        plt.tight_layout()
        #plt.show()
        plt.savefig('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/Interaction_Images/' + d1 + ',' + d2 + '.pdf', format='pdf', dpi=800)
        plt.close()
    except:
        print 'Problem'

