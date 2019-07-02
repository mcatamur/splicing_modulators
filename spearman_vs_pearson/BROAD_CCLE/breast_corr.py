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
            self.parser.add_argument('-i', '--inputFile', action = 'store', help = 'input file 1, file with ccle outlier data')
            self.parser.add_argument('-c', '--color', action = 'store', help = 'drug name', default='green')


            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input = myCommandLine.args.inputFile
        self.color = myCommandLine.args.color
        self.get_outlier_data()
        self.plot_correlation()

    def get_outlier_data(self):
        """
        Save data of outlier counts from input file
        """

        self.breast_data = {'MDA_MB_175': 4,
        'ZR_75_30': 5,
        'CAMA_1' : 8,
        'MDA_MB_134' : 13,
        'HCC_202' : 21,
        'UACC_893' : 24,
        'EFM_19' : 27,
        'SUM_190': 28,
        'EFM_129A': 42,
        'MDA_MB_361': 44,
        'HCC_1500': 45,
        'HCC_1419': 51,
        'HCC_38': 64,
        'MDA_MB_415': 64,
        'MCF_10A': 92,
        'UACC_812': 96,
        'HCC_2218': 100,
        'ZR_75_1': 110,
        'MDA_MB_435S': 115,
        '184A1': 118,
        'T47_D': 127,
        'MCF_7': 148,
        'BT_20': 177,
        'MDA_MB_435': 201,
        'BT_474': 240,
        'SK_BR_3': 300,
        'KPL_1': 327,
        'HCC_1143': 359,
        'MDA_MB_435': 201,
        'BT-474': 240,
        'SK_BR_3': 300,
        'KPL_1': 327,
        'HCC_1143': 359,
        'MDA_MB_231': 432,
        'HCC_1395': 472,
        'SUM_225': 503,
        'Hs_578T': 524,
        '184B5': 538,
        'UACC_732': 744,
        'CAL_51': 905,
        'BT_549': 1000,
        'COLO_824': 1000,
        'DU4475': 1000,
        'HCC_1187': 1000,
        'HCC_1569': 1000,
        'HCC_1806': 1000,
        'HCC_1937': 1000,
        'HCC_1954': 1000,
        'HCC_70': 1000,
        'MDA_MB_436': 1000,
        'MDA_MB_157': 1000,
        'MDA_MB_468': 1000
        }

        key_list = self.breast_data.keys()


        self.breast_plot_data = {}
        with open(self.input) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                a=row[0].split('-')
                #check if input file has overlaps:
                if a[3] in key_list:
                    self.breast_plot_data[a[3]] = int(row[1])




    def plot_correlation(self):
        """
        Now plot ic50 values vs outlier counts
        """

        #append in for loop to make sure the pairing remain in order

        sorted_dict = sorted(self.breast_data.items(), key=operator.itemgetter(1))

        x = []
        y = []
        for key in sorted_dict:
            try:
                nm_convert =  (self.breast_data[key[0]])*(10**(-9))
                log_nm_convert = math.log10(nm_convert)
                plt.scatter(log_nm_convert, self.breast_plot_data[key[0]], label=key[0],  marker = '.', s = 200, alpha = 0.6 )
                x.append(log_nm_convert)
                y.append(math.log10(self.breast_data[key[0]]))
            except KeyError:
                pass


        spearman_coefficient = spearmanr(x,y)

        pearson_coefficient = pearsonr(x,y)
        plt.text(.5,.7,'spearman: %.3f \n pval: %.3f \n pearson: %.3f \n pval: %.3f' % (spearman_coefficient[0], spearman_coefficient[1], pearson_coefficient[0], pearson_coefficient[1]), bbox={'facecolor':'w','pad':5}, \
         ha="center", va="bottom", transform=plt.gca().transAxes)




        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.07), \
          fancybox=True, shadow=True, ncol=5)
        plt.plot(x, y, self.color)
        plt.title('Palboclib_BREAST')
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
