# include electrode_distance information 
## Look at distance and number relationship
from pulse2percept.implants import ArgusII
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
def dva2ret(r_deg):
    # from pulse2percept
    """Converts visual angles (deg) into retinal distances (um)
    This function converts degrees of visual angle into a retinal distance
    from the optic axis (um) using Eq. A5 in [Watson2014]_.
    Parameters
    ----------
    r_dva : double or array-like
        Eccentricity in degrees of visual angle (dva)
    Returns
    -------
    r_um : double or array-like
        Eccentricity in microns
    """
    sign = np.sign(r_deg)
    r_deg = np.abs(r_deg)
    r_mm = 0.268 * r_deg + 3.427e-4 * r_deg ** 2 - 8.3309e-6 * r_deg ** 3
    r_um = 1e3 * r_mm
    return sign * r_um

#if it is <2, the value should be 0.
def convert_prob(row):
    if int(row['pts_number_processed'])>=2:
        return 1
    else:
        return 0


def combine_electrodes(row):
    e1 = row['pts_electrode1']
    e2 = row['pts_electrode2']
    electrodes = [e1, e2]
    electrodes.sort()
    return electrodes[0] + '-' + electrodes[1]
def get_distance(row, simulations):
    e1 =  simulations.loc[row['subject_id'],'implant'][row['pts_electrode1']]
    e2 = simulations.loc[row['subject_id'],'implant'][row['pts_electrode2']]
    #fi1nd the euclidean distance of two given electrodes
    distance = np.sqrt((e1.x -e2.x)**2 + (e1.y-e2.y)**2)
   # print('Electrode Pair %s: distance= %4.3f \n'%(row['electrode_pair'], distance))
    return distance 

def get_axonal_locations(row, model_type, simulations):
    argus = simulations.loc[row['subject_id'],'implant']
    model = simulations.loc[row['subject_id'],'axonmap']
    axons = model.grow_axon_bundles()
    e1 = argus[row['pts_electrode1']]
    e2 = argus[row['pts_electrode2']]
    a1 = model.find_closest_axon(axons, e1.x, e1.y)
    a1 = a1[0] #get the axon out of the single-item array
    a2 = model.find_closest_axon(axons, e2.x, e2.y)
    a2 = a2[0]
    
    ind_e1 = np.argmin(np.sqrt((a1[:, 0] - e1.x) ** 2 + (a1[:, 1] - e1.y) ** 2))
    ind_e2 = np.argmin(np.sqrt((a2[:, 0] - e2.x) ** 2 + (a2[:, 1] - e2.y) ** 2))

   #find the distance between e1 and a1
    if model_type == 'unidirectional': 
        # only look at the portion of the axon starting from electrode antidromically
        dist_e1a2 =  np.min(np.sqrt((a2[ind_e2:, 0] - e1.x) ** 2 + (a2[ind_e2:, 1] - e1.y) ** 2))
        ind_e1a2 =  np.argmin(np.sqrt((a2[ind_e2:, 0] - e1.x) ** 2 + (a2[ind_e2:, 1] - e1.y) ** 2))
     
        dist_e2a1 =  np.min(np.sqrt((a1[ind_e1:, 0] - e2.x) ** 2 + (a1[ind_e1:, 1] - e2.y) ** 2))
        ind_e2a1 =  np.argmin(np.sqrt((a1[ind_e1:, 0] - e2.x) ** 2 + (a1[ind_e1:, 1] - e2.y) ** 2))
        
        
        # add the sliced indexes so that you return the index of the original array
        # rather than the sliced array
        # this works and it's magic
        ind_e1a2 = ind_e1a2 + ind_e2
        ind_e2a1 = ind_e2a1 + ind_e1
    else:
        dist_e1a2 =  np.min(np.sqrt((a2[:, 0] - e1.x) ** 2 + (a2[:, 1] - e1.y) ** 2))
        ind_e1a2 =  np.argmin(np.sqrt((a2[:, 0] - e1.x) ** 2 + (a2[:, 1] - e1.y) ** 2))

        dist_e2a1 =  np.min(np.sqrt((a1[:, 0] - e2.x) ** 2 + (a1[:, 1] - e2.y) ** 2))
        ind_e2a1 =  np.argmin(np.sqrt((a1[:, 0] - e2.x) ** 2 + (a1[:, 1] - e2.y) ** 2))
    #get distance to the axon from the respective electrode
    distance_to_axon = np.minimum(dist_e1a2, dist_e2a1)

    #get distance along the axon
    #choose the small one
    if distance_to_axon == dist_e1a2:
        return {'axon':a2, 
                'axon2': a1,
                'ind_axon':ind_e1a2,
                'ind_ea' : np.argmin(np.sqrt((a2[:, 0] - e2.x) ** 2 + (a2[:, 1] - e2.y) ** 2)),      
                'ea_x' : e2.x, #electrode on the given axon
                'ea_y' : e2.y,
                'er_x' : e1.x, #other electrode (respective)
                'er_y' : e1.y,
                'phy_distance' : np.sqrt((e1.x -e2.x)**2 + (e1.y-e2.y)**2),
                'distance_to_axon': distance_to_axon}

    else:
        return {'axon' : a1,
                'axon2': a2,
                'ind_axon' : ind_e2a1,
                'ind_ea' : np.argmin(np.sqrt((a1[:, 0] - e1.x) ** 2 + (a1[:, 1] - e1.y) ** 2)),
                'ea_x' : e1.x,
                'ea_y' : e1.y, 
                'er_x' : e2.x,
                'er_y' : e2.y,
                'phy_distance' : np.sqrt((e1.x -e2.x)**2 + (e1.y-e2.y)**2),
                'distance_to_axon': distance_to_axon}

def get_axonal_distances(row, model_type, simulations):
    
    # retrieve position and distances of the row
    locations = get_axonal_locations(row, model_type, simulations)
    # returns the index of the point on the axon where electrode coincides and
    # the distance from that point to the optic nerve
    if locations['ind_ea'] <=locations['ind_axon']:
        segment_along_axon = locations['axon'][locations['ind_ea']:locations['ind_axon']]
    else:
        segment_along_axon =locations['axon'][locations['ind_axon']:locations['ind_ea']]

    # for each consecutive axong segment, we take the absolute
    # value of the difference between two points, then sum all the differences 
    # to calculate distance along the axon
    distance_along_axon = np.sum(np.abs([segment_along_axon[i + 1] - segment_along_axon[i] for i in range(len(segment_along_axon)-1)]))
    locations['distance_along_axon'] = distance_along_axon
    locations['segment_along_axon'] = segment_along_axon
    return pd.Series(locations)



def plot_distances(df, ax, flip = False):
    cmap= plt.get_cmap("Dark2") 
    ax.set_xlim([-5000,5000])
    ax.set_ylim([-5000,5000])
    #plot axon1
    ax.plot(df['axon'][:, 0], df['axon'][:, 1], c=(0.5, 0.5, 0.5))
    #plot axon 2 
    ax.plot(df['axon2'][:, 0], df['axon2'][:, 1], c=(0.5, 0.5, 0.5))


    #plot physical distance
    pplot,= ax.plot([df['ea_x'], df['er_x']], [df['ea_y'], df['er_y']], lw=3, 	color=cmap.colors[0], alpha=0.7)

    #plot axonal distance
    aplot,= ax.plot([df['er_x'], df['axon'][df['ind_axon'],0]], [ df['er_y'],  df['axon'][df['ind_axon'],1]],lw=3, color=cmap.colors[1],alpha=0.7)
    #ax.margins(tight=True)

    #plot optic nerve
    ax.add_patch(patches.Circle(df['axon'][0], radius=900, alpha=1,
                                 color='gray', zorder=10))

    #plot distance along axon
    cplot,= ax.plot(df['segment_along_axon'][:,0],df['segment_along_axon'][:,1],  lw=3, color=cmap.colors[2],alpha=0.7)

    #plot electrode1
    ax.plot( df['ea_x'], df['ea_y'],'or', markersize=10)
    #plot electrode 2
    ax.plot(df['er_x'],df['er_y'],'or', markersize=10)

    ax.legend(handles = [pplot,aplot,cplot],labels=['Physical Distance: %i' %df['phy_distance'], 
                                                    'Distance to Axon: %i'%df['distance_to_axon'], 'Distance Along Axon: %i' %df['distance_along_axon']])
    if flip:
        ax.invert.yaxis()


    
