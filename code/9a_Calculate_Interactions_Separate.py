import numpy as np
from sklearn.decomposition import PCA
from sklearn.neighbors.kde import KernelDensity

from scipy.spatial import distance as dis
from scipy.spatial.distance import mahalanobis

import matplotlib

#matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.patches import Patch

from sklearn.cluster import KMeans

import os






# Load all the drug perturbation vectors
path = '../data/Calculate_Interactions/All_Vectors_Combined.csv'
fp = open(path)

features = fp.readline().split(',')[1:]
numfeatures = len(features)

drug_perturbation_vectors = {'Batch1': {}, 'Batch2': {}}
dmso = {'Batch1': [], 'Batch2': []}

All_CLOUDs = set()
# Go threw the file and create DMSO and drug_perturbation_vectors
# DMSO per batch
# drug_perturbation_vectors per Batch with full identifier
for line in fp:
    tmp = line.strip().split(',')

    drug1, drug2 = tmp[0].split('_')[0].split('|')
    well = tmp[0].split('_')[1]
    plate = tmp[0].split('_')[2]

    values = list(np.float_(tmp[1:]))

    if int(plate) < 1315065:
        if drug1 == 'DMSO':
            dmso['Batch1'].append(values)
        drug_perturbation_vectors['Batch1'][drug1 + ',' + drug2 + ',' + plate + ',' + well] = values

    else:
        if drug1 == 'DMSO':
            dmso['Batch2'].append(values)
        drug_perturbation_vectors['Batch2'][drug1 + ',' + drug2 + ',' + plate + ',' + well] = values

        All_CLOUDs.add(drug1)

# Get list of All_CLOUDs, this number is smaller than the original number due to some drugs not being transformer correctly or other problems
All_CLOUDs.remove('DMSO')
All_CLOUDs.remove('PosCon')

All_CLOUDs = list(All_CLOUDs)
All_CLOUDs.sort()
print 'Number of drugs that have at least one correct well: %d' % len(All_CLOUDs)

# Both thresholds need to be true to set a drug as decayed during experiment
threshold_decay = 0.05
threshold_MaxDifference = 0.3

# Load all the drug decay regressions
# Created by checking the single drug responses over the different plates (there is a temporal context between plate 1 and 123)
# One is interested both in the decay as well as the maximum change e.g. if gradient between 0.1 to 0.2, still ok
# Create a dic that tells about the status of drug decay i.e. True if drug WORKED CORRECTLY
path = '../data/Calculate_Interactions/DrugDecay_Combined.csv'
fp = open(path)
fp.next()
drug_decay = {}
batch1_Failed = 0
batch2_Failed = 0
for line in fp:
    tmp = line.strip().split(',')

    batch1_decay = float(tmp[1])
    batch1_diff = float(tmp[2])

    batch2_decay = float(tmp[3])
    batch2_diff = float(tmp[4])

    batch1_Status = True
    if batch1_decay >= threshold_decay and batch1_diff >= threshold_MaxDifference:
        batch1_Status = False
        batch1_Failed += 1

    batch2_Status = True
    if batch2_decay >= threshold_decay and batch2_diff >= threshold_MaxDifference:
        batch2_Status = False
        batch2_Failed += 1

    drug_decay[tmp[0]] = {'Batch1': batch1_Status, 'Batch2': batch2_Status}
fp.close()

print 'Number of drugs that decayed in batch1: %d' % batch1_Failed
print 'Number of drugs that decayed in batch2: %d' % batch2_Failed

cutoff_min_cells = 30

well_cell_count = {}
empty_well = 0
path = '../data/Calculate_Interactions/All_CellCounts_Combined.csv'
fp = open(path)
fp.next()
for line in fp:
    tmp = line.strip().split(',')

    num_cells = float(tmp[4])

    well_cell_count[tmp[0] + ',' + tmp[1] + ',' + tmp[2] + ',' + tmp[3]] = {'Number': num_cells, 'Worked': tmp[5]}

print 'Mean number of cells in Batch1: %f' % np.mean(
    [well_cell_count[x]['Number'] for x in well_cell_count if int(x.split(',')[2]) < 1315065])
print 'Mean number of cells in Batch1: %f' % np.mean(
    [well_cell_count[x]['Number'] for x in well_cell_count if int(x.split(',')[2]) >= 1315065])

print 'Mean empty wells in Batch1: %d' % len([well_cell_count[x]['Number'] for x in well_cell_count if
                                              int(x.split(',')[2]) < 1315065 and well_cell_count[x][
                                                  'Number'] < cutoff_min_cells])
print 'Mean empty wells in Batch2: %d' % len([well_cell_count[x]['Number'] for x in well_cell_count if
                                              int(x.split(',')[2]) >= 1315065 and well_cell_count[x][
                                                  'Number'] < cutoff_min_cells])

DMSO_CellCount = {}
DMSO_CellCount['Batch1'] = np.percentile([well_cell_count[x]['Number'] for x in well_cell_count if
                                          int(x.split(',')[2]) < 1315065 and x.split(',')[0] == 'DMSO'], 90)
DMSO_CellCount['Batch2'] = np.percentile([well_cell_count[x]['Number'] for x in well_cell_count if
                                          int(x.split(',')[2]) >= 1315065 and x.split(',')[0] == 'DMSO'], 90)


def Calculate_MahalanobisDistance(k, numberTreatment):
    # Fit a PCA
    pca = PCA()
    pca.fit(k)

    # Take only as many principel components so taht 90% of the variance can be explained
    variances = pca.explained_variance_ratio_
    total = 0
    use_components = 0
    for v in variances:
        total = total + v
        use_components += 1
        if total > 0.9:
            break

    if use_components <= 1:
        use_components = 2



    pca = PCA(n_components=use_components)
    pca.fit(k)

    transformedX = pca.transform(k)

    weightedPCA = np.multiply(transformedX, pca.explained_variance_ratio_)

    pca1 = list(weightedPCA[:, 0])
    pca2 = list(weightedPCA[:, 1])

    treatments = weightedPCA[0:numberTreatment]
    control = weightedPCA[numberTreatment:]

    # Contains the main treatment vector
    x = []
    for i in range(0, use_components):
        x.append(np.mean(treatments[:, i]))

    # Contains the mean untreated vector
    u = []
    for i in range(0, use_components):
        u.append(np.mean(control[:, i]))

    S = np.cov(weightedPCA, rowvar=0)

    x = np.array(x)
    u = np.array(u)

    distance_calc = mahalanobis(x, u, np.linalg.inv(S))

    # original = weightedPCA.copy()
    treatment_row = weightedPCA[0].copy()

    all_mahala_distances = []
    for i in range(0, 1000):
        np.random.shuffle(weightedPCA)

        # if np.array_equal(original, weightedPCA):
        if np.array_equal(treatment_row, weightedPCA[0]):
            continue

        treatments = weightedPCA[0:numberTreatment]
        control = weightedPCA[numberTreatment:]

        x = []
        for i in range(0, use_components):
            x.append(np.mean(treatments[:, i]))

        u = []
        for i in range(0, use_components):
            u.append(np.mean(control[:, i]))




        covariance_treatment = np.cov(treatments, rowvar=0)
        covariance_control = np.cov(control, rowvar=0)

        weighted_covariance_treatment = covariance_treatment * (
        float(len(treatments)) / (len(treatments) + len(control)))
        weighted_covariance_control = covariance_control * (float(len(control)) / (len(treatments) + len(control)))

        S = weighted_covariance_treatment + weighted_covariance_control

        S = np.cov(weightedPCA, rowvar=0)

        x = np.array(x)
        u = np.array(u)
        distance_rand = mahalanobis(x, u, np.linalg.inv(S))
        all_mahala_distances.append(distance_rand)

    # print len([x for x in all_mahala_distances if x >= distance_calc]) / float(len(all_mahala_distances))

    # caluclate empirical pvalue
    mp = len([x for x in all_mahala_distances if x >= distance_calc]) / float(len(all_mahala_distances))
    if mp > 1:
        mp = 1

    return distance_calc, mp, pca1, pca2



cloud = 'CLOUD001'

for b in ['Batch1', 'Batch2']:

    # Create a list with Treatments (drug vectors) and DMSO
    k = [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if
         cloud + ',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]
    numberTreatment = len(k)
    k.extend(dmso[b])




    distance_calc, mp, pca1, pca2 = Calculate_MahalanobisDistance(k, numberTreatment)

    if mp < 0.01:
        color = ['#b2182b'] * numberTreatment
    else:
        color = ['#2166ac'] * numberTreatment

    for i in range(1, len(k)):
        color.append('grey')

    plt.scatter(pca1[numberTreatment:], pca2[numberTreatment:], c=color[numberTreatment:], alpha=0.4)
    plt.scatter(pca1[0:numberTreatment], pca2[0:numberTreatment], c=color[0:numberTreatment], alpha=0.4)
    plt.legend(['Mahalanobis Distance: %.2f\n mp = %.2e' % (distance_calc, mp)])
    plt.show()

    plt.close()





exit()
Single_Drug_Significance = {}

if os.path.isfile('../results/Calculate_Interactions/Singles/Overview.csv'):
    fp = open('../results/Calculate_Interactions/Singles/Overview.csv', 'r')
    fp.next()
    for line in fp:
        tmp = line.strip().split(',')
        Single_Drug_Significance[tmp[0]] = {
            'Batch1': {'Mahalanobis_Distance': float(tmp[1]), 'MP_Value': float(tmp[2])},
            'Batch2': {'Mahalanobis_Distance': float(tmp[3]), 'MP_Value': float(tmp[4])}}
else:

    for cloud in All_CLOUDs:
        print cloud
        Single_Drug_Significance[cloud] = {}
        for b in ['Batch1', 'Batch2']:

            # Create a list with Treatments (drug vectors) and DMSO
            k = [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if
                 cloud + ',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]
            numberTreatment = len(k)
            k.extend(dmso[b])

            distance_calc, mp, pca1, pca2 = Calculate_MahalanobisDistance(k, numberTreatment)

            if mp < 0.01:
                color = ['#b2182b'] * numberTreatment
            else:
                color = ['#2166ac'] * numberTreatment

            for i in range(1, len(k)):
                color.append('grey')

            Single_Drug_Significance[cloud][b] = {'Mahalanobis_Distance': distance_calc, 'MP_Value': mp}

            plt.scatter(pca1[numberTreatment:], pca2[numberTreatment:], c=color[numberTreatment:], alpha=0.4)
            plt.scatter(pca1[0:numberTreatment], pca2[0:numberTreatment], c=color[0:numberTreatment], alpha=0.4)
            plt.legend(['Mahalanobis Distance: %.2f\n mp = %.2e' % (distance_calc, mp)])
            # plt.show()
            plt.savefig('../results/Calculate_Interactions/Singles/' + cloud + '_' + b + '.png')
            plt.close()

    fp_out = open('../results/Calculate_Interactions/Singles/Overview.csv', 'w')
    fp_out.write('Drug,Batch1_Mahalanobis_Distance,Batch1_MP_Value,Batch2_Mahalanobis_Distance,Batch2_MP_Value\n')
    for cloud in All_CLOUDs:
        fp_out.write(cloud + ',' + str(Single_Drug_Significance[cloud]['Batch1']['Mahalanobis_Distance']) + ',' + str(
            Single_Drug_Significance[cloud]['Batch1']['MP_Value']) + ',' + str(
            Single_Drug_Significance[cloud]['Batch2']['Mahalanobis_Distance']) + ',' + str(
            Single_Drug_Significance[cloud]['Batch2']['MP_Value']) + '\n')
    fp_out.close()

# print Single_Drug_Significance

print 'Number of significant drugs in Batch1: %d' % len(
    [x for x in Single_Drug_Significance if Single_Drug_Significance[x]['Batch1']['MP_Value'] < 0.001])
print 'Number of significant drugs in Batch2: %d' % len(
    [x for x in Single_Drug_Significance if Single_Drug_Significance[x]['Batch2']['MP_Value'] < 0.001])

# check if drugs have enough replicates
enough_replicates = {}
for cloud in All_CLOUDs:
    enough_replicates[cloud] = {}
    for b in ['Batch1', 'Batch2']:

        # Create a list with Treatments (drug vectors) and DMSO
        k = [x for x in drug_perturbation_vectors[b] if cloud + ',DMSO' in x]

        too_few_cells = 0
        for n in k:
            tmp = n.split(',')

            if well_cell_count[tmp[0] + ',' + tmp[1] + ',' + tmp[2] + ',' + tmp[3]]['Number'] < cutoff_min_cells:
                too_few_cells += 1

        if (len(k) - too_few_cells) < 3:
            enough_replicates[cloud][b] = False
        else:
            enough_replicates[cloud][b] = True

# print enough_replicates
print 'Number of drugs with too few replicates in batch1: %d' % len(
    [x for x in enough_replicates if enough_replicates[x]['Batch1'] == False])
print 'Number of drugs with too few replicates in batch2: %d' % len(
    [x for x in enough_replicates if enough_replicates[x]['Batch2'] == False])

# list of drugs that kill too many cells to do morphological analysis with them
killing_drugs = {'Batch1': [], 'Batch2': []}

# this drugs must be removed from further analysis
drugs_to_remove = {'Batch1': [], 'Batch2': []}

for cloud in All_CLOUDs:
    for b in ['Batch1', 'Batch2']:

        enough_replicates_status = enough_replicates[cloud][b]
        drug_decay_status = drug_decay[cloud][b]

        if enough_replicates_status and drug_decay_status:
            continue
        else:
            # print 'Bad Drug (%s): %s    Decay?: %s' %(b,cloud, str(not drug_decay_status))

            if drug_decay_status:
                number_replicates = len([x for x in drug_perturbation_vectors[b] if cloud + ',DMSO' in x])

                if number_replicates >= 3:
                    # print '==> can be used a "killing drug"'
                    killing_drugs[b].append(cloud)
                else:
                    drugs_to_remove[b].append(cloud)
            else:
                drugs_to_remove[b].append(cloud)

print 'Drugs to remove:'
print drugs_to_remove

print 'Cytotoxic drugs/non usable morphological features:'
print killing_drugs


# Some Easy Outlier detection
def reject_outliers_2(data, m=6.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / (mdev if mdev else 1.)
    return [data[i] for i in range(0, len(data)) if s[i] < m]


# Actual Math for calculating the DDIs
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            angle_between((1, 0, 0), (1, 0, 0))
            0.0
            angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def calculate_vector_math_v2(a, b, c):
    '''
    calculate the amount of single a, single b, and the 'surprise factor)
    :param a: vector a (single)
    :param b: vector b (single)
    :param c: vector c (combination)
    :return: alpha, beta and gamma (part of vector a, b and suprise)
    '''

    if sum(a) != 0 and sum(b) != 0:

        if angle_between(a, b) <= 0.5:
            # h = c/(a+b)

            A = np.array([[np.dot(a, a), np.dot(a, b)], [np.dot(a, b), np.dot(b, b)]])
            h = np.array([np.dot(a, c), np.dot(b, c)])

            alpha = h[0] / (A[0][0] + A[1][0])
            beta = h[1] / (A[0][1] + A[1][1])

            n = alpha * a + beta * b - c

            gamma = np.linalg.norm(n)
            return str(alpha), str(beta), str(gamma)

        elif angle_between(a, b) == 3.141592653589793:
            A = np.array([[np.dot(a, a), np.dot(a, b)], [np.dot(a, b), np.dot(b, b)]])
            h = np.array([np.dot(a, c), np.dot(b, c)])

            alpha = h[0] / (A[0][0] + abs(A[1][0]))
            beta = h[1] / (abs(A[0][1]) + A[1][1])

            n = alpha * a + beta * b - c

            gamma = np.linalg.norm(n)

            return str(alpha), str(beta), str(gamma)

    try:

        if sum(c) != 0 and sum(a) == 0 and sum(b) == 0:
            return '1', '1', str(dis.euclidean([0] * len(c), c))
        elif sum(c) == 0 and sum(a) == 0 and sum(b) == 0:
            return '1', '1', '0'
        elif sum(c) == 0 and sum(a) == 0:
            return '1', '0', '0'
        elif sum(c) == 0 and sum(b) == 0:
            return '0', '1', '0'

        else:
            # Matrix equation
            A = np.array([[np.dot(a, a), np.dot(a, b)], [np.dot(a, b), np.dot(b, b)]])

            h = np.array([np.dot(a, c), np.dot(b, c)])

            if A[0][0] == 0 and h[0] == 0:  # one vector zero, so combination can be only 1dim

                beta = h[1] / A[1][1]
                n = beta * b - c
                gamma = np.linalg.norm(n)

                if (len(list(b)) - list(b).count(0) > 2) or np.linalg.norm(b) > 0.5:
                    return '1.0', str(beta), str(gamma)
                else:
                    return '1', '1', str(dis.euclidean([0] * len(c), c))
            elif A[1][1] == 0 and h[1] == 0:
                alpha = h[0] / A[0][0]
                n = alpha * a - c
                gamma = np.linalg.norm(n)

                if len(list(a)) - list(a).count(0) > 2 or np.linalg.norm(a) > 0.5:
                    return str(alpha), '1', str(gamma)
                else:
                    return '1', '1', str(dis.euclidean([0] * len(c), c))
            elif h[0] == 0 and h[1] == 0:
                gamma = np.linalg.norm(c)
                return '0.0', '0.0', str(gamma)
            elif A[0][0] != 0 and h[0] == 0 and h[1] == 0:
                gamma = np.linalg.norm(c)
                return '0.0', '1.0', str(gamma)
            elif A[1][1] != 0 and h[0] == 0 and h[1] == 0:
                gamma = np.linalg.norm(c)
                return '1.0', '0', str(gamma)

            p = np.linalg.solve(A, h)
            # orthogonal vector
            n = p[0] * a + p[1] * b - c

            distance = np.linalg.norm(n)
            # check
            # print('dot product of a and c: %.4f' %(np.dot(a,n)))
            # print('dot product of b and c: %.4f' %(np.dot(b,n)))
            # print('distance: %.3f' %(distance))
            return str(p[0]), str(p[1]), str(distance)
    except:
        return 'Error', 'Error', 'Error'


def Create_Single_Drug_Vectors():
    drug_vectors = {'Batch1': {}, 'Batch2': {}}
    for cloud in All_CLOUDs:
        for b in ['Batch1', 'Batch2']:
            if Single_Drug_Significance[cloud][b]['MP_Value'] < 0.01:

                # change this func
                drug_wells = [x for x in drug_perturbation_vectors[b] if
                              cloud + ',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]

                drug_well_values = []
                for well in drug_wells:
                    drug_well_values.append(drug_perturbation_vectors[b][well])

                # s_Vectors =  [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if cloud+',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]

                keep_names = []
                keep_values = []
                for s1_a, label_a in zip(drug_well_values, drug_wells):
                    tmp = []
                    for s1_b, label_b in zip(drug_well_values, drug_wells):
                        # if s1_a != s1_b:
                        sim = 1 - dis.cosine(s1_a, s1_b)
                        tmp.append(sim)
                    if np.mean(tmp) > 0.3:
                        keep_names.append(label_a)
                        keep_values.append(s1_a)
                        # similarity.append(tmp)

                drug_vectors[b][cloud] = {}
                for val, lab in zip(keep_values, keep_names):
                    drug_vectors[b][cloud][lab] = val

    return drug_vectors


drug_vectors_WithoutOutliers = Create_Single_Drug_Vectors()


alreadyExisting = [f.split('.')[0] for f in os.listdir('../results/Calculate_Interactions/MC_Scores/') if os.path.isfile(os.path.join('../results/Calculate_Interactions/MC_Scores/', f))]



numbers = range(1,len(All_CLOUDs))



for my_Number in numbers:


    cloud1 = All_CLOUDs[my_Number]
    if cloud1 in alreadyExisting:
        continue



    #fp_out = open('../results/Calculate_Interactions/Bliss_Scores/' + cloud1 + '.csv', 'w')
    #fp_out.write('Drug1,Drug2,Batch,Plate,Well,S1,S2,C,Expected,Difference\n')

    #fp_out2 = open('../results/Calculate_Interactions/MC_Scores/' + cloud1 + '.csv', 'w')
    #fp_out2.write('Drug1,Drug2,Batch,Plate,Well,S1_Mahalanobis,S1_MPValue,S1_Norm,S2_Mahalanobis,S2_MPValue,S2_Norm,Combi_Mahalanobis,Combi_MPValue,Combi_Norm_Norm,Alpha,Beta,Gamma,Alpha_Zero,Beta_Zero,Gamma_Zero,Mahalanobis_Distance_To_NI,MPValue_To_NI,Mahalanobis_Distance_To_NI_NoNoise,MPValue_To_NI_NoNoise\n')

    Combination_Drug_Significance = {}

    print cloud1
    for cloud2 in All_CLOUDs:
        if cloud1 > cloud2:

            cloud1 = 'CLOUD192'
            cloud2 = 'CLOUD065'

            combination_well = [x for x in well_cell_count if cloud1 + ',' + cloud2 in x or cloud2 + ',' + cloud1 in x]

            if len(combination_well) == 0:
                print 'Problem!!!'

            plate = combination_well[0].split(',')[2]
            well = combination_well[0].split(',')[3]

            if int(plate) < 1315065:
                b = 'Batch1'
            else:
                b = 'Batch2'

            if cloud1 in drugs_to_remove[b] or cloud2 in drugs_to_remove[b]:
                #fp_out2.write(cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + ','.join(
                #    ['DrugProblem'] * 19) + '\n')
                continue

            single1_replicate_wells = [x for x in drug_perturbation_vectors[b] if cloud1 + ',DMSO' in x]
            single2_replicate_wells = [x for x in drug_perturbation_vectors[b] if cloud2 + ',DMSO' in x]

            combi_Number = well_cell_count[combination_well[0]]['Number']
            combi_Worked = well_cell_count[combination_well[0]]['Worked']

            s1_cellCount = [well_cell_count[x]['Number'] for x in single1_replicate_wells if
                            well_cell_count[x]['Worked'] == 'TRUE']
            s2_cellCount = [well_cell_count[x]['Number'] for x in single2_replicate_wells if
                            well_cell_count[x]['Worked'] == 'TRUE']

            if combi_Worked == 'FALSE':
                #fp_out2.write(cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + ','.join(
                #    ['TransferProblem'] * 19) + '\n')
                continue

            else:

                s1_cellCount = reject_outliers_2(s1_cellCount, m=4)
                s2_cellCount = reject_outliers_2(s2_cellCount, m=4)

                s1_cytotoxicity = (DMSO_CellCount[b] - np.mean(s1_cellCount)) / DMSO_CellCount[b]
                if s1_cytotoxicity < 0:
                    s1_cytotoxicity = 0
                s2_cytotoxicity = (DMSO_CellCount[b] - np.mean(s2_cellCount)) / DMSO_CellCount[b]
                if s2_cytotoxicity:
                    s2_cytotoxicity = 0

                comb_cytotoxicity = (DMSO_CellCount[b] - combi_Number) / DMSO_CellCount[b]
                if comb_cytotoxicity < 0:
                    comb_cytotoxicity = 0

                Bliss_expected = s1_cytotoxicity + s2_cytotoxicity - s1_cytotoxicity * s2_cytotoxicity
                if Bliss_expected > 1:
                    Bliss_expected = 1

                Bliss_Difference = comb_cytotoxicity - Bliss_expected

                #fp_out.write(cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + str(
                #    s1_cytotoxicity) + ',' + str(s2_cytotoxicity) + ',' + str(comb_cytotoxicity) + ',' + str(
                #    Bliss_expected) + ',' + str(Bliss_Difference) + '\n')
                # if combi_Number < cutoff_min_cells:



                if cloud1 in killing_drugs[b] or cloud2 in killing_drugs[b]:
                    #fp_out2.write(cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + ','.join(
                    #    ['CytotoxicDrug'] * 19) + '\n')
                    continue



                    #################
                    # FROM HERE actual morphology
                    # 1.) check if the two singles and/or combination are significant morphological changers
                    # 2.) depending on that calculate interactions
                    ###########


                    # Check if combination is significant
                if drug_perturbation_vectors[b].has_key(combination_well[0]) == False:
                    #fp_out2.write(cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + ','.join(
                    #    ['ProblemWithCombinationWell'] * 19) + '\n')
                    continue
                combi_Vector = drug_perturbation_vectors[b][combination_well[0]]

                k = [combi_Vector]
                numberTreatment = len(k)
                k.extend(dmso[b])

                combination_distance_calc, combination_mp, _, _ = Calculate_MahalanobisDistance(k, numberTreatment)
                Combi_VectorNorm = np.linalg.norm(combi_Vector)

                s1_Mahalanobis = Single_Drug_Significance[cloud1][b]['Mahalanobis_Distance']
                s1_MValue = Single_Drug_Significance[cloud1][b]['MP_Value']

                s2_Mahalanobis = Single_Drug_Significance[cloud2][b]['Mahalanobis_Distance']
                s2_MValue = Single_Drug_Significance[cloud2][b]['MP_Value']

                no_interaction = False
                if (combination_mp >= 0.05 or combination_distance_calc < 6) and (
                        s1_MValue >= 0.01 or s1_Mahalanobis < 6) and (s2_MValue >= 0.01 or s2_Mahalanobis < 6):
                    alpha = 1
                    beta = 1
                    gamma = 0
                    alpha_0 = 1
                    beta_0 = 1
                    gamma_0 = 0

                    combination_distance_calc_toNI = 0
                    combination_mp_toNI = 1

                    combination_distance_calc_toNI_noNoise = 0
                    combination_mp_toNI_noNoise = 1

                    no_interaction = True

                s1_significant = False
                if s1_MValue < 0.01 and s1_Mahalanobis > 6:
                    s1_Vectors = drug_vectors_WithoutOutliers[b][cloud1].values()

                    s1_significant = True
                else:
                    s1_Vectors = [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if
                                  cloud1 + ',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]

                val = np.array(s1_Vectors)

                s1_MeanVector = []
                for i in range(0, numfeatures):
                    s1_MeanVector.append(np.mean(val[:, i]))

                s1_Mean_VectorNorm = np.linalg.norm(s1_MeanVector)

                s2_significant = False
                if s2_MValue < 0.01 and s2_Mahalanobis > 6:
                    s2_Vectors = drug_vectors_WithoutOutliers[b][cloud2].values()
                    s2_significant = True
                else:
                    s2_Vectors = [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if
                                  cloud2 + ',DMSO' in x and well_cell_count[x]['Number'] > cutoff_min_cells]

                val = np.array(s2_Vectors)

                s2_MeanVector = []
                for i in range(0, numfeatures):
                    s2_MeanVector.append(np.mean(val[:, i]))
                s2_Mean_VectorNorm = np.linalg.norm(s2_MeanVector)

                comb_significant = False
                if combination_mp < 0.05 and combination_distance_calc > 6:
                    comb_significant = True

                if no_interaction == False:

                    s1_MeanVector_toCalculate = np.array(s1_MeanVector)
                    s2_MeanVector_toCalculate = np.array(s2_MeanVector)
                    comb_MeanVector_toCalculate = np.array(combi_Vector)

                    alpha, beta, gamma = calculate_vector_math_v2(s1_MeanVector_toCalculate, s2_MeanVector_toCalculate,
                                                                  comb_MeanVector_toCalculate)

                    if s1_significant == True:
                        s1_MeanVector_toCalculate = np.array(s1_MeanVector)
                    else:
                        s1_MeanVector_toCalculate = np.array([0] * numfeatures)

                    if s2_significant == True:
                        s2_MeanVector_toCalculate = np.array(s2_MeanVector)
                    else:
                        s2_MeanVector_toCalculate = np.array([0] * numfeatures)

                    if comb_significant == True:
                        comb_MeanVector_toCalculate = np.array(combi_Vector)
                    else:
                        comb_MeanVector_toCalculate = np.array([0] * numfeatures)

                    alpha_0, beta_0, gamma_0 = calculate_vector_math_v2(s1_MeanVector_toCalculate,
                                                                        s2_MeanVector_toCalculate,
                                                                        comb_MeanVector_toCalculate)

                    mu, sigma = 0, 0.05
                    numberTreatment = 1

                    possible_results = []

                    # Create all possible vector sums (theoretical points of NI)
                    first_vector = []
                    if len(s1_Vectors) >= len(s2_Vectors):
                        first_vector = s1_Vectors
                        second_Vector = s2_Vectors
                        num_Clusters = len(s1_Vectors)
                    else:
                        first_vector = s2_Vectors
                        second_Vector = s1_Vectors
                        num_Clusters = len(s2_Vectors)

                    for s1 in first_vector:
                        tmp = []
                        for s2 in second_Vector:
                            possible_results.append(np.array(s1) + np.array(s2))

                    k = list(possible_results)
                    k.insert(0, combi_Vector)
                    combination_distance_calc_toNI_noNoise, combination_mp_toNI_noNoise, _, _ = Calculate_MahalanobisDistance(
                        k, numberTreatment)

                    noise_vectors = []
                    for i in range(1, 10):
                        noise = np.random.normal(mu, sigma, [len(possible_results), numfeatures])
                        added_noise_points = possible_results + noise
                        noise_vectors.extend(added_noise_points)

                    possible_results.extend(noise_vectors)

                    kmeans = KMeans(n_clusters=num_Clusters, n_init=100, init='random')
                    kmeans.fit(possible_results)
                    y_kmeans = kmeans.predict(possible_results)
                    centers = set(y_kmeans)

                    results_MPVal_Mahala = []

                    for c in centers:
                        cluster_values = []
                        for v, l in zip(possible_results, y_kmeans):
                            if l == c:
                                cluster_values.append(v)

                        k = cluster_values
                        k.insert(0, combi_Vector)
                        combination_distance_calc_toNI, combination_mp_toNI, _, _ = Calculate_MahalanobisDistance(k,
                                                                                                                  numberTreatment)

                        results_MPVal_Mahala.append((combination_distance_calc_toNI, combination_mp_toNI))

                        # plt.scatter(pca1[1:], pca2[1:], c='grey', alpha=0.4)
                        # plt.scatter(pca1[0:1], pca2[0:1], c='red', alpha=0.4)
                        # plt.show()

                    combination_mp_toNI = max([x[1] for x in results_MPVal_Mahala])

                    combination_distance_calc_toNI = min(
                        [x[0] for x in results_MPVal_Mahala if x[1] == combination_mp_toNI])

                    # Check the mahalanobis distance between the combination and the theoretical NIs

                    k = list(possible_results)
                    k.insert(0, combi_Vector)

                    # print len(possible_results)

                    # Make a PCA plot if it is a significant interaction
                    if (combination_mp_toNI < 0.05 and combination_distance_calc_toNI > 3) and (
                            combination_mp_toNI_noNoise < 0.05 and combination_distance_calc_toNI_noNoise > 3):

                        number_vectorsums = len(k)

                        color = ['#b2182b'] * numberTreatment
                        for i in range(1, len(k)):
                            color.append('blue')

                        DMSO_Wells = [drug_perturbation_vectors[b][x] for x in drug_perturbation_vectors[b] if
                                      'DMSO,None,' + str(plate) in x and well_cell_count[x]['Number'] > cutoff_min_cells]
                        k.extend(DMSO_Wells)

                        for i in range(0, len(DMSO_Wells)):
                            color.append('grey')

                        k.extend(s1_Vectors)

                        for i in range(0, len(s1_Vectors)):
                            color.append('#e08214')

                        k.extend(s2_Vectors)

                        for i in range(0, len(s2_Vectors)):
                            color.append('#5aae61')

                        pca = PCA(n_components=2)
                        pca_vals = pca.fit_transform(k)

                        if s1_significant:
                            cloud1_name = cloud1 + '(+)'
                        else:
                            cloud1_name = cloud1 + '(-)'

                        if s2_significant:
                            cloud2_name = cloud2 + '(+)'
                        else:
                            cloud2_name = cloud2 + '(-)'

                        legend_elements = [Patch(facecolor='#b2182b', lw=1,
                                                 label='Combination\nMahalanobis Distance: %.2f\n mp = %.2e' % (
                                                 combination_distance_calc_toNI, combination_mp_toNI)),
                                           Patch(facecolor='blue',
                                                 label='Vector_Sums'),
                                           Patch(facecolor='grey',
                                                 label='DMSO'),
                                           Patch(facecolor='#e08214',
                                                 label=cloud1_name),
                                           Patch(facecolor='#5aae61',
                                                 label=cloud2_name)]

                        plt.scatter(pca_vals[:, 0][1:number_vectorsums], pca_vals[:, 1][1:number_vectorsums],
                                    c=color[1:number_vectorsums], alpha=0.4)

                        plt.scatter(pca_vals[:, 0][number_vectorsums:], pca_vals[:, 1][number_vectorsums:],
                                    c=color[number_vectorsums:], alpha=0.6)
                        plt.scatter(pca_vals[:, 0][0], pca_vals[:, 1][0], c=color[0], alpha=1)
                        # Single_Drug_Significance[cloud][b] = {'Mahalanobis_Distance':distance_calc,'MP_Value':mp}
                        # plt.legend(['Combination\nMahalanobis Distance: %.2f\n mp = %.2e' % (combination_distance_calc_toNI, combination_mp_toNI)])

                        plt.legend(handles=legend_elements, loc='best', prop={'size': 7})

                        # plt.scatter(pca1[1:], pca2[1:], c=color[1:], alpha=0.4)
                        # plt.scatter(pca1[0:1], pca2[0:1], c=color[0:1], alpha=0.4)
                        # plt.legend(['DMSO',cloud1_name,cloud2_name,'Vector_Sums','Combination\nMahalanobis Distance: %.2f\n mp = %.2e' % (combination_distance_calc_toNI,combination_mp_toNI)])
                        # plt.show()
                        plt.savefig('../results/Calculate_Interactions/Combinations/' + cloud1 + '_' + cloud2 + '.png')
                        plt.close()
                '''
                fp_out2.write(
                    cloud1 + ',' + cloud2 + ',' + str(b) + ',' + plate + ',' + well + ',' + str(s1_Mahalanobis) + ',' + str(
                        s1_MValue) + ',' + str(s1_Mean_VectorNorm) + ',' + str(s2_Mahalanobis) + ',' + str(
                        s2_MValue) + ',' + str(s2_Mean_VectorNorm) + ',' + str(combination_distance_calc) + ',' + str(
                        combination_mp) + ',' + str(Combi_VectorNorm) + ',' + str(alpha) + ',' + str(beta) + ',' + str(
                        gamma) + ',' + str(alpha_0) + ',' + str(beta_0) + ',' + str(gamma_0) + ',' + str(
                        combination_distance_calc_toNI) + ',' + str(combination_mp_toNI) + ',' + str(
                        combination_distance_calc_toNI_noNoise) + ',' + str(combination_mp_toNI_noNoise) + '\n')
                
                '''

    #fp_out.close()
    #fp_out2.close()

    print 'Created all interactions'

print 'Make Summary File'
print 'MC_Scores:'
MC_Score_Files = [f for f in os.listdir('../results/Calculate_Interactions/MC_Scores/') if os.path.isfile(os.path.join('../results/Calculate_Interactions/MC_Scores/', f))]
MC_Score_Files.sort()
fp_out = open('../results/Calculate_Interactions/All_MC_Scores.csv','w')
fp_out.write('Drug1,Drug2,Batch,Plate,Well,S1_Mahalanobis,S1_MPValue,S1_Norm,S2_Mahalanobis,S2_MPValue,S2_Norm,Combi_Mahalanobis,Combi_MPValue,Combi_Norm_Norm,Alpha,Beta,Gamma,Alpha_Zero,Beta_Zero,Gamma_Zero,Mahalanobis_Distance_To_NI,MPValue_To_NI,Mahalanobis_Distance_To_NI_NoNoise,MPValue_To_NI_NoNoise\n')
for file_name in MC_Score_Files:
    fp = open('../results/Calculate_Interactions/MC_Scores/'+file_name,'r')
    fp.next()
    for line in fp:
        fp_out.write(line)
    fp.close()
fp_out.close()
print 'Done'

print 'Bliss_Scores:'
Bliss_Score_Files = [f for f in os.listdir('../results/Calculate_Interactions/Bliss_Scores/') if os.path.isfile(os.path.join('../results/Calculate_Interactions/Bliss_Scores/', f))]
Bliss_Score_Files.sort()
fp_out = open('../results/Calculate_Interactions/All_Bliss_Scores.csv','w')
fp_out.write('Drug1,Drug2,Batch,Plate,Well,S1,S2,C,Expected,Difference\n')
for file_name in Bliss_Score_Files:
    fp = open('../results/Calculate_Interactions/Bliss_Scores/'+file_name,'r')
    fp.next()
    for line in fp:
        fp_out.write(line)
    fp.close()
fp_out.close()
print 'Done'