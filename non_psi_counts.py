""" author: mcatamur@ucsc.edu """


"""
Sample Command: python IQR_tissue_outlier.py -i juncBase_table.txt/tsv
"""

#!/usr/bin/env python3


import csv, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
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
            self.parser.add_argument('-t','--tissue_type', type=str, action = 'store', help = 'tissue_type to analyze', default='LUNG')
            self.parser.add_argument('-dPSI','--delta_PSI', type=int, action = 'store', help = 'delta PSI threshold', default=10.0)

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.tissue = myCommandLine.args.tissue_type
        self.dPSI = myCommandLine.args.delta_PSI
        self.findTissueIndex(self.inputFile)
        #self.buildData_tissue(self.inputFile)
        #self.plotEventsBar()
        #self.plotEventsBar_allTissue()
        #self.plotEventsBar_oneTissue()


    def findTissueIndex(self, inputFile):
        self.nonpsi_ccle = {}
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for ccle in range(11,len(header)-1): #index of cell lines
                self.nonpsi_ccle[header[ccle]] = 0

        with open(self.inputFile) as file:
            with open(self.inputFile) as file:
                tsvReader = csv.reader(file, dialect = 'excel-tab')
                header = next(tsvReader, None)
                for row in tsvReader:
                    for ccle in range(11,len(header)-1): #index of cell lines
                        try:
                            if row[ccle] == 'NA':
                                self.nonpsi_ccle[header[ccle]] += 1
                        except ValueError:
                            pass

        sorted_dict = sorted(self.nonpsi_ccle.items(), key=operator.itemgetter(1))
        for key in sorted_dict:
            print(key[0], key[1])


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
