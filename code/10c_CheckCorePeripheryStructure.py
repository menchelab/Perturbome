'''
This part because of cpalgorithm only existing for python 3.x.x was written in PYTHON3
Take care to use the right version for this script!
'''

import cpalgorithm as cp
import networkx as nx
import matplotlib
from matplotlib import pylab as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import numpy as np

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'



DPI = nx.read_gml('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/Networks/DPI_Network_Complete.gml')
print('Number of edges: %d' %len(DPI.edges()))
print('Number of nodes: %d' %len(DPI.nodes()))

#Models: BE, MINRES, SBM, Rossa, Surprise

model = 'MINRES'

#only needed for continues models like Rossa
cut = 0.15

algorithm = cp.MINRES()
algorithm.detect(DPI)

#dictionary of nodes, telling to which core-periphery pair the node belongs (if only single core
#-periphery structure then then each key has the value 0)
c = algorithm.get_pair_id()
#Dictionary with node as keys telling if the node is a core node (1) or periphery node(0)
x = algorithm.get_coreness()

# this function calculates the significance that the network is a core perihpery network by comparing it to
# 1000 random networks (configuration network)
sig_c, sig_x, significant, p_values = cp.qstest(c, x, DPI, algorithm, num_of_rand_net=500)
#sig_c same as c but only for significant assiciations (else None value)
#sig_x same as x but only for significant assiciations (else None value)
#boolean True/False if network is core perihpery structures
#p_values for the individual core-periphery structures
print(sig_x)


print('Significant Core-Periphery Structure?')
print(significant)
print('Pvalues: ')
print(p_values)


if model != 'Rossa':
    print('Core Nodes: %d' %len([node for node in sig_x if sig_x[node] == 1]))
    print([node for node in sig_x if sig_x[node] == 1])


    print('Periphery Nodes: %d' %len([node for node in sig_x if sig_x[node] == 0]))
    print([node for node in sig_x if sig_x[node] == 0])

else:
    print('here')
    print('Core Nodes: %d' % len([node for node in sig_x if sig_x[node] > cut]))
    print([node for node in sig_x if sig_x[node]> cut])

    print('Periphery Nodes: %d' % len([node for node in sig_x if sig_x[node]  <= cut]))
    print([node for node in sig_x if sig_x[node] <= cut])

fp_out = open('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/CorePerihpery/'+model+'/Result.csv','w')
fp_out.write('Significant,' + str(significant[0]) +'\n')
fp_out.write('PValue,%.4f\n' %p_values[0])
if model != 'Rossa':
    fp_out.write('CoreNodes,%s\n' %';'.join([node for node in sig_x if sig_x[node] == 1]))
    fp_out.write('PerihperyNodes,%s' %';'.join([node for node in sig_x if sig_x[node] == 0]))
else:
    fp_out.write('CoreNodes,%s\n' % ';'.join([node for node in sig_x if sig_x[node] > cut]))
    fp_out.write('PerihperyNodes,%s' % ';'.join([node for node in sig_x if sig_x[node] <= cut]))
fp_out.close()

degree = nx.degree(DPI)

top10 = np.percentile([x[1] for x in degree],90)


core_degrees = []
perihpery_degrees = []
if model != 'Rossa':

    for tuple in degree:
        if sig_x[tuple[0]] == 1:
            core_degrees.append(tuple[1])
        elif sig_x[tuple[0]] == 0:
            perihpery_degrees.append(tuple[1])
else:
    for tuple in degree:
        if sig_x[tuple[0]] > cut:
            core_degrees.append(tuple[1])
        elif sig_x[tuple[0]] <= cut:
            perihpery_degrees.append(tuple[1])

plt.hist(core_degrees, color = '#FA9900', density=True)
plt.hist(perihpery_degrees, color = '#A1999D', density=True)

# Create the formatter using the function to_percent. This multiplies all the
# default labels by 100, making them all percentages
formatter = FuncFormatter(to_percent)
# Set the formatter
plt.gca().yaxis.set_major_formatter(formatter)
plt.axvline(top10, ls='--', color='grey')
plt.legend(['10% Highest Degree'])
#plt.show()
plt.savefig('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/CorePerihpery/'+model+'/Histogram.pdf')
plt.close()




bplot = sns.boxplot(data=[perihpery_degrees,core_degrees], orient='h',showfliers=True)
interaction_types = ['Perihpery','Core']
interaction_colors = ['#A1999D','#FA9900']
color_dict = dict(zip(interaction_types, interaction_colors))
for i in range(0, 2):
    mybox = bplot.artists[i]
    mybox.set_facecolor(color_dict[interaction_types[i]])
plt.axvline(top10, ls='--', color='grey')
plt.legend(['10% Highest Degree'])
#plt.show()
plt.savefig('../results/Create_DPI_Network/DPI_iS3_pS7_abMAD2_gP100/CorePerihpery/'+model+'/Boxplot.pdf')
plt.close()