""" author: mcatamur@ucsc.edu """


"""
Get ic50 values for the cell lines in the input file and plot them against
outlier counts.
"""

#!/usr/bin/env python3

import csv, time, argparse
start_time = time.time()
import operator
import math
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests




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
            self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'output file')

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input1 = myCommandLine.args.inputFile1
        self.input2 = myCommandLine.args.inputFile2
        self.outputFile = myCommandLine.args.outputFile
        self.all_data = {}
        self.get_outlier_data()
        self.get_ic50()
        self.multipleTestCorrection()


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
                self.ccle_list[row[0]] = [row[2], int(row[3]), int(row[1])]

            #print(self.ccle_list)

        #row[0] -  cell line name
        #row[1] - outlier counts
        #row[2] - tissue type
        #row[3] - index

    def get_ic50(self):
        """
        parse file for ic50 values
        """
        drug_list = []
        with open(self.input2, encoding='ISO-8859-1') as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                drug_list.append(row[1])

            #print (drug_list)

        for drug in drug_list:
            self.getSpearmanr(drug)

        #print(self.all_data)

    def getSpearmanr(self, drug):
        with open(self.input2, encoding='ISO-8859-1') as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            lymph_data = []
            cns_data = []
            lung_data = []
            for row in tsvReader:
                if row[1] == drug:
                    for ccle in self.ccle_list:
                        if  self.ccle_list[ccle][0] == 'LUNG':
                            try:
                                lung_data.append([math.log10((float(row[self.ccle_list[ccle][1]]))), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                            except ValueError:
                                try:
                                    a = row[self.ccle_list[ccle][1]]
                                    b = float(a[1:])
                                    if b == 20:
                                        lung_data.append([math.log10(20.01), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    elif b == 0.001221:
                                        lung_data.append([math.log10(0.001220), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    else:
                                        lung_data.append([math.log10(b), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                except ValueError:
                                    pass


                        ##################################################################################################

                        if  self.ccle_list[ccle][0] == 'CENTRAL_NERVOUS_SYSTEM':
                            try:
                                cns_data.append([math.log10((float(row[self.ccle_list[ccle][1]]))), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                            except ValueError:
                                try:
                                    a = row[self.ccle_list[ccle][1]]
                                    b = float(a[1:])
                                    if b == 20:
                                        cns_data.append([math.log10(20.01), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    elif b == 0.001221:
                                        cns_data.append([math.log10(0.001220), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    else:
                                        cns_data.append([math.log10(b), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                except ValueError:
                                    pass


                        ##################################################################################################

                        if  self.ccle_list[ccle][0] == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE':
                            try:
                                lymph_data.append([math.log10((float(row[self.ccle_list[ccle][1]]))), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                            except ValueError:
                                try:
                                    a = row[self.ccle_list[ccle][1]]
                                    b = float(a[1:])
                                    if b == 20:
                                        lymph_data.append([math.log10(20.01), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    elif b == 0.001221:
                                        lymph_data.append([math.log10(0.001220), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                    else:
                                        lymph_data.append([math.log10(b), self.ccle_list[ccle][2], self.ccle_list[ccle][0]])
                                except ValueError:
                                    pass

                    else:
                        pass

        self.build_data(drug, lung_data)
        self.build_data(drug, lymph_data)
        self.build_data(drug, cns_data)


    def build_data(self, drug, data):
        key_name = drug + '_' + data[0][2]

        x=[]
        y=[]

        for item in data:
            x.append(item[0])
            y.append(item[1])


        spearman_coefficient = spearmanr(x,y)
        #print(spearman_coefficient[1])
        #pearson_coefficient = pearsonr(x,y)


        if math.isnan(spearman_coefficient[1]):
            pass
        else:
            self.all_data[key_name] = [float(spearman_coefficient[1]), float(spearman_coefficient[0])]

        self.p_values = {}
        for key in self.all_data:
            self.p_values[key] = self.all_data[key][0]



    def multipleTestCorrection(self):
        sorted_dict1 = sorted(self.all_data.items(), key=operator.itemgetter(0))
        p_values =[]

        sorted_dict2 = sorted(self.p_values.items(), key=operator.itemgetter(0))

        for key in sorted_dict2:
            p_values.append(key[1])

        #a = fdr_correction(p_values, alpha=0.05, method='indep')
        a, b, c, d = multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        print(multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False))
        with open(self.outputFile, mode ='w') as out_file:
            tsvWriter = csv.writer(out_file, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            tsvWriter.writerow(['drug : tissue_type' + '\t' + 'spearman p-value' + '\t' + 'spearman corr val' + '\t' + 'adjusted p-value'])
            for i in range(0, len(b)):
                #print (sorted_dict1[i])

                tsvWriter.writerow([sorted_dict1[i][0] + '\t' + str(sorted_dict1[i][1][0]) + '\t' + str(sorted_dict1[i][1][1]) + '\t' +  str(b[i])])

















def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
