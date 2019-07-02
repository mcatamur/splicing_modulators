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
            self.parser.add_argument('-i1', '--inputFile1', action = 'store', help = 'input file 1, file with ccle outlier data')
            self.parser.add_argument('-i2', '--inputFile2', action = 'store', help = 'input file 2, file with ic50 data')
            self.parser.add_argument('-t', '--tissueType', action = 'store', help = 'tissueType', default='HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
            self.parser.add_argument('-d', '--drug', action = 'store', help = 'drug name', default='PHA-665752')
            self.parser.add_argument('-c', '--color', action = 'store', help = 'drug name', default='green')



            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input1 = myCommandLine.args.inputFile1
        self.input2 = myCommandLine.args.inputFile2
        self.tissue = myCommandLine.args.tissueType
        self.drug = myCommandLine.args.drug
        self.color = myCommandLine.args.color
        self.get_outlier_data()
        self.get_ic50()
        self.plot_correlation()

    def get_outlier_data(self):
        """
        Save data of outlier counts from input file
        """
        self.ccle_list = {}
        self.outlier_data = {}
        with open(self.input1) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                self.ccle_list[row[0]] = [row[2], int(row[3])]
                self.outlier_data[row[0]] = int(row[1])

        #row[0] -  cell line name
        #row[1] - outlier counts
        #row[2] - tissue type
        #row[3] - index

    def get_ic50(self):
        """
        parse file for ic50 values
        """
        self.ic50 = {}
        with open(self.input2, encoding='ISO-8859-1') as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                if row[1] == self.drug:
                    for ccle in self.ccle_list:
                        if  self.ccle_list[ccle][0] == self.tissue:
                            try:
                                self.ic50[ccle] = math.log10((float(row[self.ccle_list[ccle][1]])))
                            except ValueError:
                                try:
                                    a = row[self.ccle_list[ccle][1]]
                                    b = float(a[1:])
                                    if b == 20:
                                        self.ic50[ccle]= math.log10(20.01)
                                    elif b == 0.001221:
                                        self.ic50[ccle] = math.log10(0.001220)
                                    else:
                                        self.ic50[ccle]= math.log10(b)
                                        print(math.log10(b))
                                except ValueError:
                                    pass

                else:
                    pass


    def plot_correlation(self):
        """
        Now plot ic50 values vs outlier counts
        """

        #append in for loop to make sure the pairing remain in order
        outlier = {}

        for item in self.ic50:
            outlier[item] = self.outlier_data[item]

        #first sort the dictionary
        sorted_dict = sorted(self.ic50.items(), key=operator.itemgetter(1))

        x = []
        y = []
        for item in sorted_dict:
            x.append(item[1])
            y.append(self.outlier_data[item[0]])
            plt.scatter(item[1], self.outlier_data[item[0]], label=item[0],  marker = '.', s = 200, alpha = 0.6 )

        spearman_coefficient = spearmanr(x,y)

        pearson_coefficient = pearsonr(x,y)
        plt.text(.5,.7,'spearman: %.3f \n pval: %.3f \n pearson: %.3f \n pval: %.3f' % (spearman_coefficient[0], spearman_coefficient[1], pearson_coefficient[0], pearson_coefficient[1]), bbox={'facecolor':'w','pad':5}, \
         ha="center", va="bottom", transform=plt.gca().transAxes)




        #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.07), \
        #  fancybox=True, shadow=True, ncol=5)
        plt.plot(x, y, self.color)
        plt.title('%s: %s' % (self.drug, self.tissue))
        plt.xlabel('log(IC50)')
        plt.ylabel('outlier counts')
        plt.savefig('%s_%s' % (self.drug, self.tissue), dpi=1200)
        plt.show()






def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
