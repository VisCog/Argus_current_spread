import pulse2percept as pulse2percept
from scipy.spatial.distance import cdist
from pulse2percept.models import ScoreboardModel, AxonMapModel
from pulse2percept.implants import ArgusII
import pandas as pd

def get_sims(subjectdata):
	tsample = 0.01
	img_thresh = 60
	c_amp = 330


	axon_models = []
	score_models = []
	subjects = subjectdata.index

	implants = []
	for subject in subjects:
	        x_center= subjectdata.loc[subject, 'implant_x']
	        y_center = subjectdata.loc[subject, "implant_y"]
	        rot = subjectdata.loc[subject, 'implant_rot']
	        eye = subjectdata.loc[subject, 'eye']
	        implant = subjectdata.loc[subject, 'implant_type_str']
	        loc_od_x = subjectdata.loc[subject, 'loc_od_x']
	        loc_od_y = subjectdata.loc[subject, 'loc_od_y']
	        rho = subjectdata.loc[subject, 'rho']
	        no_axons = 500
	        axlambda = subjectdata.loc[subject, 'lambda']

	        smodel= ScoreboardModel(rho=rho,
	                 xrange=(-20, 20), xystep=0.25,
	                 yrange=(-15, 15))
	        score_models.append(smodel.build())


	        amodel= AxonMapModel(rho=rho, axlambda=axlambda, ignore_pickle = True,
	                 eye= eye, loc_od=(loc_od_x,loc_od_y),n_axons=no_axons,
	                 n_jobs=-1, scheduler='threading',
	                 thresh_percept=0, verbose=True,
	                 xrange=(-20, 20), xystep=0.25,
	                 yrange=(-15, 15))
	        axon_models.append(amodel.build())


	        #generate model
	        argus = ArgusII(x=x_center, y=y_center,rot=rot,eye=eye)
	        implants.append(argus)
	        
	simulations = pd.DataFrame()

	simulations['scoreboard']= score_models
	simulations['axonmap']= axon_models
	simulations['implant']=implants
	simulations['subjects'] = subjects
	simulations = simulations.set_index('subjects')

	return simulations









