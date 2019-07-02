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
import matplotlib.pyplot as plt
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
            self.parser.add_argument('-i', '--inputFile', action = 'store', help = 'input file, file with ccle outlier data from step 1')
            #self.parser.add_argument('-i2', '--inputFile2', action = 'store', help = 'input file 2, file with ic50 data')
            #self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'output file')
            self.parser.add_argument('-d', '--drug', action = 'store', help = 'compound name')

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input = myCommandLine.args.inputFile
        #self.input2 = myCommandLine.args.inputFile2
        #self.outputFile = myCommandLine.args.outputFile
        self.drug = myCommandLine.args.drug
        self.all_data = {}
        self.get_outlier_data()
        self.checkTissue()
        self.multipleTestCorrection()

        """
        if self.drug == 'Docetaxel':
            self.plot_correlation('DOCETAXEL_LUNG')
            self.plot_correlation('DOCETAXEL_CENTRAL_NERVOUS_SYSTEM')
            self.plot_correlation('DOCETAXEL_BONE_AND_SOFT_TISSUE')
        elif self.drug == 'Staurosporine':
            self.plot_correlation('Staurosporine_LUNG')
        """

    def get_outlier_data(self):
        """
        Save data of outlier counts from input file
        """
        self.ccle_list = {}
        #row_num = 0
        with open(self.input) as file:

            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                #row_num += 1
                self.ccle_list[row[0]] = [int(row[1]), (row[2]), float(row[3])]

                #print (row_num)

            #print(self.ccle_list)

        #row[0] -  cell line name
        #row[1] - outlier counts
        #row[2] - tissue type
        #row[3] - ic50
        #print(self.all_data)
    def checkTissue(self):
        overlap = {}

        for key in self.ccle_list:
            try:
                overlap[self.ccle_list[key][1]] += 1
            except KeyError:
                overlap[self.ccle_list[key][1]] = 1



        for key in overlap:
            if overlap[key] >= 10:
                self.getSpearmanr(key)


    def getSpearmanr(self, key):
        key_data = []

        for ccle in self.ccle_list:
            if  self.ccle_list[ccle][1] == key:
                    key_data.append([math.log10(self.ccle_list[ccle][2]), self.ccle_list[ccle][0], self.ccle_list[ccle][1], ccle])
            else:
                pass


        self.build_data(key_data)

        if self.drug == 'Docetaxel':
            if key == 'BONE_AND_SOFT_TISSUE':
                self.plot_correlation(key_data)
            elif key == 'LUNG':
                self.plot_correlation(key_data)
            elif key == 'CENTRAL_NERVOUS_SYSTEM':
                self.plot_correlation(key_data)
        elif self.drug == 'Staurosporine':
            if key == 'LUNG':
                self.plot_correlation(key_data)
        elif self.drug == 'Imatinib':
            if key == 'DLBC':
                self.plot_correlation(key_data)





    def build_data(self, data):
        #print (data)
        #print (data[0][2])
        key_name = self.drug + '_' + data[0][2]

        x=[]
        y=[]

        for item in data:
            x.append(item[0])
            y.append(item[1])


        spearman_coefficient = spearmanr(x,y)
        #print(spearman_coefficient[0])
        #print(spearman_coefficient[1])
        pearson_coefficient = pearsonr(x,y)
        #print(pearson_coefficient[1])
        #print(pearson_coefficient[1])



        if math.isnan(spearman_coefficient[1]):
            pass
        else:
            self.all_data[key_name] = [float(spearman_coefficient[1]), float(spearman_coefficient[0]), float(pearson_coefficient[1]), float(pearson_coefficient[0])]


        self.p_values = {}
        self.p_values_pearson = {}

        for key in self.all_data:
            self.p_values[key] = self.all_data[key][0]
            self.p_values_pearson[key] = self.all_data[key][2]

    def plot_correlation(self, data):
        """
        Now plot ic50 values vs outlier counts
        """
        #append in for loop to make sure the pairing remain in order

        key_name = self.drug + '_' + data[0][2]
        sort = sorted(data, key=lambda x:x[0])
        #print(sort)


        x = []
        y = []
        for item in sort:
            x.append(item[0])
            y.append(item[1])
            plt.scatter(item[0], item[1], label=item[3],  marker = '.', s = 200, alpha = 0.6 )

        spearman_coefficient = spearmanr(x,y)

        pearson_coefficient = pearsonr(x,y)
        plt.text(.5,.7,'spearman: %.3f \n pval: %.3f \n pearson: %.3f \n pval: %.3f' % (spearman_coefficient[0], spearman_coefficient[1], pearson_coefficient[0], pearson_coefficient[1]), bbox={'facecolor':'w','pad':5}, \
         ha="center", va="bottom", transform=plt.gca().transAxes)




        #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.07), \
        #  fancybox=True, shadow=True, ncol=5)
        plt.plot(x, y, 'green')
        plt.title(key_name)
        plt.xlabel('log(IC50)')
        plt.ylabel('outlier counts')
        plt.savefig(key_name, dpi=600)
        plt.show()

    def multipleTestCorrection(self):
        #print(self.all_data)

        output = self.drug + '_' + 'GDSC_2.txt'

        sorted_dict1 = sorted(self.all_data.items(), key=operator.itemgetter(0))

        p_values = []
        p_values2 = []

        sorted_dict2 = sorted(self.p_values.items(), key=operator.itemgetter(0))
        sorted_dict3 = sorted(self.p_values_pearson.items(), key=operator.itemgetter(0))
        for key in sorted_dict2:
            p_values.append(key[1])
        for key in sorted_dict3:
            p_values2.append(key[1])


        #print(sorted_dict2)
        #a = fdr_correction(p_values, alpha=0.05, method='indep')
        a, b, c, d = multipletests(p_values, alpha=0.25, method='fdr_bh', is_sorted=False, returnsorted=False)

        e, f, g, h = multipletests(p_values2, alpha=0.25, method='fdr_bh', is_sorted=False, returnsorted=False)

        #print(multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False))
        with open(output, mode ='w') as out_file:
            #tsvWriter = csv.writer(out_file, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            out_file.write('drug:tissue_type' + '\t' + 'spearman p-value' + '\t' + 'spearman corr val' + '\t' + 'adjusted spearman p-value' + '\t' + 'pearson p-value' + '\t' + 'pearson corr val' + '\t' + 'adjusted pearson p-value' + '\n')
            for i in range(0, len(b)):
                #print (sorted_dict1[i])
                out_file.write(sorted_dict1[i][0] + '\t' + str(sorted_dict1[i][1][0]) + '\t' + str(sorted_dict1[i][1][1]) + '\t' +  str(b[i]) + '\t' + str(sorted_dict1[i][1][2]) + '\t' + str(sorted_dict1[i][1][3]) + '\t' + str(f[i]) + '\n')


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
