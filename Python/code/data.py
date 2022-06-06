import os
import pandas as pd
import os.path as op
import logging

class processing:
    def __init__(self):
        #set up logger
        logging.getLogger().handlers = []
        self.logger = logging.getLogger(__name__)
        self.logger.addHandler(logging.StreamHandler())

        # Set up root logger for debug file
        formatstr = '%(asctime)s  %(message)s'
        logging.basicConfig(level=logging.INFO,
                            format=formatstr,
                            filename='.img_selection.log',
                            filemode='w')
    def get_two_point(self, outpath, scriptpath, filename):


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
    # fixes the problem with string integer conflict for PTS number
    # modified from the initial preprocessing code to only account for double
    # experiments
    def process_ptsnumber_in_double(self, row, change23=True):
        #disable Jupyter Notebook handlers
        
        val = row['pts_number']

        #if value is a string
        if isinstance(val, str):
            if self._is_digit(val): #a float or int 
                if change23:
                    if int(val) >2:
                        return 2
                    else:
                        return int(float(val))
                else:
                        return int(float(val))
            else: #when val is string and not a digit
                #check some common cases 
                if val =='nothing' or val =='-':
                    return 0
                else:
                    user_input = input('There is a string instead of a number,\n'
                                       'Please fix %s\n'
                                       'Description is: %s'%(val, row['pts_description']))
                    inp = self._add_pts_number(user_input)
                    print(inp)
                    return inp
        #if value is nan
        elif pd.isnull(val):
            user_input = input('No number is added for file %s'
                               '\n, \n could you fix this problem by looking at description \n %s?'
                               %(row['pts_file'], row['pts_description'] ))

            inp = self._add_pts_number(user_input)
            print(inp)
            return inp

            #if value is float
        elif isinstance(val, float) or isinstance(val, int):
            if change23:
                if int(val) >2:
                    return 2
                else:
                    return int(val)
            else:
                return int(val)
        else:
            print('This is weird')

    #helper function
    def _is_digit(self, x):
        try:
            float(x)
            return True
        except ValueError:
            return False

    #helper function for user input data for pts_number_processing
    def _add_pts_number(self, user_input):
        if pd.isnull(user_input):
                self.logger.info('OK, you are not adding any file, this may cause a problem later.')
                return np.nan
        elif user_input.isdigit():
            user_input = int(user_input)
            #check if it is within meaningful limits
            if user_input >=0 or user_input<=100:
                self.logger.info('You are adding the integer  %i' %user_input)
                return user_input

            else:
                self.logger.info('Your value does not make sense.')
                self.logger.info('No number is added')
                return np.nan
        else:
            logger.info('Not a valid type, nothing is added.')
            return np.nan
