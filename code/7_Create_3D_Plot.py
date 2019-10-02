#USE MAYAVI TO CREATE A SHINY 3D PLOT ILLUSTRATING THE CLOUD PERTURBATION SPACE
#

#import all necessary python libraries
import networkx as nx
import numpy as np
from mayavi import mlab


def hex_to_rgb(value):
    '''
    Easy function that converts a HEX color code into an RGB color code.
    Returns the corosponding RGB color code.
    '''
    value = value.lstrip('#')
    lv = len(value)

    return tuple(int(value[i:i+lv/3], 16)/255.0 for i in range(0, lv, lv/3))


# ===============================================================
def plot_3D(posfile='crys_coords_small_200.dat', colfile='None', draw_edges=0, point_scale=0.25):


    # Define the drug module color. Colorize:
    # C03(Diuretics), 
    # P01(Anti0Protozoal), 
    # N04(Anti Parkinson),
    # N06(Psychoanaleptics) and 
    # M04(Antigout) as example ATC classes; 
    # other ATC classes  stay grey
    colors = {}
    if colfile != 'None':
        fp = open(colfile)
        fp.next()
        for line in fp:
            tmp = line.strip().split(',')

            if 'C03' in tmp[3]:
                colors[tmp[0]] =  hex_to_rgb('#B0FFFF')
            elif 'P01' in tmp[3]:
                colors[tmp[0]] =  hex_to_rgb('#7EC96E')
            elif 'N04' in tmp[3]:
                colors[tmp[0]] =  hex_to_rgb('#257BD2')
            elif 'N06' in tmp[3]:
                colors[tmp[0]] =  hex_to_rgb('#B2CAEA')
            elif 'M04' in tmp[3]:
                colors[tmp[0]] =  hex_to_rgb('#F9BE69')
            else:
                colors[tmp[0]] = hex_to_rgb('#bdbdbd')


    # open the position file and readlines
    f = open(posfile, 'r')
    lines = f.readlines()

    #lists containing the x,y,z,size and color property for each drug module
    l_x = []
    l_y = []
    l_z = []
    l_s = []
    l_c = []
        
    # go through all drug modules
    for line in lines[1:]:
    
        #extract CLOUD name
        drug = line.strip().split(',')[0]

        #extract position and size
        x = float(line.strip().split(',')[1])
        y = float(line.strip().split(',')[2])
        z = float(line.strip().split(',')[3])
        s = float(line.strip().split(',')[4])

        #add to lists
        l_x.append(x)
        l_y.append(y)
        l_z.append(z)
        l_s.append(s)

        #get ATC color from previously create dictionary
        l_c.append(colors[drug])
        
    #close file pointer
    f.close()

    #Go through all points and plot each point individually with the correspoinding attributes
    for x, y, z, s,c in zip(l_x, l_y, l_z, l_s,l_c):
        pts = mlab.points3d(x, y, z,
                            scale_factor=s,
                            colormap='jet',
                            color=c,  # change (with list)
                            resolution=200,  # 20 fast but not as pretty, 200 pretty but not fast
                            opacity=1)

    mlab.show()
    #mlab.savefig('test.pdf')



if __name__ == '__main__':
    # --------------------------------------------------------
    #
    # Create the drug module bubble plot
    #
    # --------------------------------------------------------

    # Get the x,y,z posititions for the drug modules
    fn_xyz = '../results/ATC_Analysis/Bubble_Positions_TargetsOnly.csv'
    
    # Find the ATC associations for each CLOUD drug
    colorFile = '../data/ATC_Analysis/CLOUD_to_ATC.csv'
    
    # Create the 3D bubble plot
    plot_3D(posfile=fn_xyz, colfile=colorFile, draw_edges=0, point_scale=1)
    sys.exit(0)



