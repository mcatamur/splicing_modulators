""" author: mcatamur@ucsc.edu """


"""
Produce a list of cell lines we have data for.
"""

#!/usr/bin/env python3

import csv, time, argparse
start_time = time.time()
import pandas as pd

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
            self.parser.add_argument('-i1', '--inputFile1', action = 'store', help = 'input file 1, file with outlier data')
            self.parser.add_argument('-i2', '--inputFile2', action = 'store', help = 'input file 2, file with ic50 data')
            #self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'output file')
            self.parser.add_argument('-d', '--drug', action = 'store', help = 'compund_name')


            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input1 = myCommandLine.args.inputFile1
        self.input2 = myCommandLine.args.inputFile2
        self.drug = myCommandLine.args.drug
        #self.output = myCommandLine.args.outputFile
        self.buildData_file1()
        self.buildData_file2()
        self.outputFile()

    def buildData_file1(self):
        """
        looks through file with outlier counts and saves it in a dictionary.
        """
        self.raw_ccle_list = {}
        with open(self.input1) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                a = row[0].split("-")
                self.raw_ccle_list[a[3]] = [row[1], row[2]]
        #print(self.raw_ccle_list)


    def buildData_file2(self):
        """
        looks through file with ic50 values and if an outlier count value is available
        for that cell line then it is appended to a final list of cell lines
        """
        self.final_ccle_list = {}

        keys = list(self.raw_ccle_list.keys())

        #overlap = []
        #tissue_count = {}
        #row_num =0
        with open(self.input2, encoding='ISO-8859-1') as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                #row_num += 1
                a = row[0].replace("-", "_")
                if a in keys:
                    self.final_ccle_list[a] = [self.raw_ccle_list[a][0], self.raw_ccle_list[a][1], row[4]]
                #print(row_num)
        """
        tissue_type = {}
        for item in self.final_ccle_list:
            try:
                tissue_type[self.final_ccle_list[item][1]] += 1
            except KeyError:
                tissue_type[self.final_ccle_list[item][1]] = 1

        print (len(self.final_ccle_list.keys()))
        print (tissue_type)
        """

    def outputFile(self):
        output = self.drug + '_' + 'GDSC_out_step1_2.txt'
        with open(output, mode ='w') as out:
            out.write('ccle'+ '\t' + 'outlier' + '\t' + 'tissue' + '\t' + 'ic50' + '\n')
            for key in self.final_ccle_list:
                out.write(key + '\t' + str(self.final_ccle_list[key][0]) + '\t' + self.final_ccle_list[key][1] + '\t' + str(self.final_ccle_list[key][2]) + '\n')




def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
