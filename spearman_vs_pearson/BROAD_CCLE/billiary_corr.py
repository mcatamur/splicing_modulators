""" author: mcatamur@ucsc.edu """


"""
Get ic50 values for the cell lines in the input file and plot them against
outlier counts.
"""

#!/usr/bin/env python3

import csv, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
import operator
import math
import numpy as np
from matplotlib.font_manager import FontProperties
from scipy.stats import spearmanr
from scipy.stats import pearsonr




class CommandLine() :
    '''
    Takes input file and options.
    '''
    def __init__(self, inOpts=None) :
            '''
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            '''
            self.parser = argparse.ArgumentParser(description = 'outputAnalyzer.py arguments',
                                                 epilog = 'For more help contact: ',
                                                 add_help = True, #default is True
                                                 prefix_chars = '-',
                                                 usage = '%(prog)s [options] -option1[default]'
                                                 )
            self.parser.add_argument('-c', '--color', action = 'store', help = 'drug name', default='green')


            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.color = myCommandLine.args.color
        self.get_outlier_data()
        self.plot_correlation()

    def get_outlier_data(self):
        """
        Save data of outlier counts from input file
        """

        self.billiary_data_ic50 = {
        'Hep_3B': 3.49,
        'JHH_7': 0.30,
        'SNU_387': 0.16,
        'SNU_449': 0.13,
        'Li_7': 0.047,
        'JHH_2': 0.022,
        'SNU_475': 0.009,
        'SNU_423': 0.009,
        'SNU_398': 0.009,
        'huH_1': 0.009,
        'HuH_7': 0.009,
        'PLC_PRF_5': 0.009,
        'SK_HEP_1': 0.009
        }

        key_list = self.billiary_data_ic50.keys()

        self.billiary_data_outlier = {
        'Hep_3B': 8791,
        'JHH_7': 8496,
        'SNU_387': 7837,
        'SNU_449': 8079,
        'Li_7': 6941,
        'JHH_2': 8539,
        'SNU_475': 8025,
        'SNU_423': 7775,
        'SNU_398': 9879,
        'huH_1': 11561,
        'HuH_7': 7752,
        'PLC_PRF_5': 8957,
        'SK_HEP_1': 7893
        }



    def plot_correlation(self):
        """
        Now plot ic50 values vs outlier counts
        """

        #append in for loop to make sure the pairing remain in order

        sorted_dict = sorted(self.billiary_data_ic50.items(), key=operator.itemgetter(1))

        x = []
        y = []
        for key in sorted_dict:
            try:
                print(self.billiary_data_ic50[key[0]])
                print((self.billiary_data_ic50[key[0]])/(10**(6)))
                micro_convert =  (self.billiary_data_ic50[key[0]])/(10**(6))
                log_micro_convert = math.log10(micro_convert)
                plt.scatter(log_micro_convert, self.billiary_data_outlier[key[0]], label=key[0],  marker = '.', s = 200, alpha = 0.6 )
                x.append(log_micro_convert)
                y.append(self.billiary_data_outlier[key[0]])
            except KeyError:
                pass


        spearman_coefficient = spearmanr(x,y)

        pearson_coefficient = pearsonr(x,y)
        plt.text(.5,.7,'spearman: %.3f \n pval: %.3f \n pearson: %.3f \n pval: %.3f' % (spearman_coefficient[0], spearman_coefficient[1], pearson_coefficient[0], pearson_coefficient[1]), bbox={'facecolor':'w','pad':5}, \
         ha="center", va="bottom", transform=plt.gca().transAxes)

        #plt.autoscale(enable=True, axis='both', tight=None)
        plt.xlim(9e-09, 2e-06)




        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.07), \
          fancybox=True, shadow=True, ncol=5)
        plt.plot(x, y, self.color)
        plt.title('Palboclib_BILLIARY')
        plt.xlabel('log(IC50)')
        plt.ylabel('outlier counts')
        plt.savefig('Palboclib_BREAST')
        plt.show()






def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
