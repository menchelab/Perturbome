import pandas
import MySQLdb
#from PIL import Image
from matplotlib import pylab as plt
import numpy as np
#from scipy import misc
#from skimage import exposure
import os

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
#from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest

def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


db = MySQLdb.connect("menchelabdb.int.cemm.at","root","cqsr4h","ImageAnalysisDDI" )

# DEFINE BATCH TO ANALYSE:
batch_ = 1


 # Get features
string = "select COLUMN_NAME from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME='DPN1018Batch"+str(batch_)+"Per_Image'"
all_features = list(pandas.read_sql(string, con=db)['COLUMN_NAME'])

i = 0
#Remove MetadataFeatures such as height of image
features = []
for f in all_features:

    if  'ImageQuality' in f and ('Max' not in f):
        features.append(f)
        print str(i)+'\t'+f +'\tQuality'
        i+=1
    if 'Image' not in f:
        print str(i) + '\t' + f + '\tMeasurement'
        i += 1
features.sort()


string = 'select ImageNumber,' +','.join(features)+' from DPN1018Batch'+str(batch_)+'Per_Image;'

ImageQuality = pandas.read_sql(string, con=db)



x = ImageQuality.loc[:, features].values
y = ImageQuality.loc[:,['ImageNumber']].values

x = StandardScaler().fit_transform(x)



print 'Train Random Forest'
clf = IsolationForest(n_jobs=1,n_estimators=500, max_samples='auto',behaviour='new', contamination='auto',random_state=1000)
clf.fit(x)
y_pred_train = clf.predict(x)
outliers = y_pred_train
print 'Finished Training'


vectorlength = np.linalg.norm(x,axis=1)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)

principalDf = pandas.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])


finalDf = pandas.concat([principalDf, ImageQuality[['ImageNumber']]], axis = 1)


print pca.explained_variance_ratio_
print pca.explained_variance_
exit()
names = finalDf['ImageNumber'].values


mean_length = np.mean(vectorlength)
standard_deviation = np.std(vectorlength)
#plt.hist(vectorlength,bins='auto')
#plt.show()
#plt.close()

colors = []
ok_image = 0
bad_image = 0
bad_images = []
bad_images_values = []
for length,o,n in zip(vectorlength,outliers,names):
#for length,n in zip(vectorlength,names):

    if o == -1:
        colors.append('red')
        bad_image += 1
        bad_images.append(n)
        bad_images_values.append(length)
    else:
        colors.append('grey')
        ok_image += 1


print 'Number of OK images: %d' %ok_image
print 'Number of bad images: %d' %bad_image

ensure_dir('../results/DPN1018Batch'+str(batch_)+'Per_Image/BadImages/BadImages.csv')
# create file containing all measurements and thresholds
fp_out = open('../results/DPN1018Batch'+str(batch_)+'Per_Image/BadImages/BadImages.csv', 'w')
for b,v in zip(bad_images,bad_images_values):
    fp_out.write(str(b)+','+str(v)+'\n')
fp_out.close()
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('ImageQuality PCA', fontsize = 20)

ax.scatter(finalDf.loc[:,'principal component 1'], finalDf.loc[:,'principal component 2'],c = colors, alpha=0.1)


ax.grid()

plt.show()


