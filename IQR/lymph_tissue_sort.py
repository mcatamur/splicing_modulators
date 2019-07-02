"""
Created on Fri Jun 23 22:54:28 2017

@author: Carmelle
"""
import csv, time, argparse
start_time = time.time()
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
import numpy as np
import math
from collections import defaultdict
import operator

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
            self.parser.add_argument('-i', '--inputFile', action = 'store', help = 'input file')
            self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'first output file')
#            self.parser.add_argument('-p', '--pValue', type=float, action = 'store', help='desired p value', default=0.01)
#            self.parser.add_argument('-sE', '--specificEvent', action = 'store', nargs='?', const=True, default=False, help='Event specific scatter plot')
#            self.parser.add_argument('-e','--eventType', type=str, action = 'store', help = 'Select an alternative splicing event type. \n Omit space; replace with underscore.', default = 'cassette')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class binaryOutliers() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.outputFile = myCommandLine.args.outputFile
        self.CCLE_dict = {}
        self.CCLE_gene_dict = {}
        self.buildData(self.inputFile)






    def buildData(self, inputFile):
        """
        Build dictionaries in this method for stacked bar plot and volcano scatter plot
        """
        lymph = {}
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader)
            for row in tsvReader:
                if row[2] == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE':
                    lymph[row[0]] = int(row[1])


        #sorted keys
        print(lymph)
        mean = sum(lymph.values())/len(lymph)



        with open(self.outputFile, mode ='w') as out_file1:
            tsvWriter = csv.writer(out_file1, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            sorted_dict = sorted(lymph.items(), key=operator.itemgetter(1))
            tsvWriter.writerow(['mean:', mean])
            tsvWriter.writerow(['median:', sorted_dict[len(sorted_dict)//2]])
            for key in sorted_dict:
                tsvWriter.writerow([key[0], key[1]])




def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = binaryOutliers(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
