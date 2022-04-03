import os
import pandas as pd
import os.path as op
def get_two_point(outpath, scriptpath, filename):


	#subject specific information
	subjectdata = pd.read_csv(op.join(outpath,'subjects-ALL.csv'),index_col='subject_id')
	n_subjects = len(subjectdata)


	editlevel = 'drawings_filled'
	subjectdata['scale'] = pd.Series([0.5] * n_subjects, index=subjectdata.index)
	subjectdata= subjectdata.drop(['28-102', '28-103', '28-106', '28-107', '28-109'])
	subjects = subjectdata.index



	data = pd.read_csv(op.join(outpath,filename))

	# We are interested in the conditions where we have patient data for all
	#select only double x2thresh experiments
	doubles = ['drawing_double_2Thr_1']
	double_x2 = data[data['name'].isin(doubles)] # double electrode at double threshold


	#calculate false positives for each patient
	#catch trials
	catch = double_x2[double_x2['pts_amp1']==0]
	catch = catch.reset_index()

	#remove catch trials from further analysis
	double_x2 = double_x2[double_x2['pts_amp1']!=0]
	double_x2.shape

	
	#remove nans from pts_number_processed
	#this is only for 12-005 where the last two trials were cut off
	double_x2= double_x2[double_x2['pts_number_processed'].notna()]
	double_x2 = double_x2[double_x2['error'].isna()]

	#remove the irrelevant error column
	double_x2 = double_x2.drop(['error'], axis =1)
	double_x2.shape


	to_drop = []
	def unique_cols(df):
	    a = df.to_numpy() # df.values (pandas<0.24)
	    return (a[0] == a).all(0)
	cols = unique_cols(double_x2)
	cols
	len(cols[cols])

	double_x2.columns[cols]

	to_drop = double_x2.columns[cols].values
	double_x2 = double_x2.drop(to_drop, axis =1)

	return double_x2, catch, subjectdata
