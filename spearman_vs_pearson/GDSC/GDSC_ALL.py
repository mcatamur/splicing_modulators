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
import os
plt.style.use('BME163')





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
            self.parser.add_argument('-i1', '--inputFile1', action = 'store', help = 'input file, file with ccle outlier data from step 1')
            self.parser.add_argument('-i2', '--inputFile2', action = 'store', help = 'input file 2, file with ic50 data')
            self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'output file')
            #self.parser.add_argument('-d', '--drug', action = 'store', help = 'compound name')

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.input1 = myCommandLine.args.inputFile1
        self.input2 = myCommandLine.args.inputFile2
        self.output = myCommandLine.args.outputFile

        os.system("touch %s" % (self.output))

        with open(self.output, 'w') as self.out_file:
            self.out_file.write('drug:tissue_type' + '\t' + 'spearman corr val' +
                        '\t' + 'spearman p-val' + '\t' + 'adjusted spearman p-val' +
                        '\t' + 'pearson corr val' + '\t' + 'pearson p-val' + '\t' +
                        'adjusted pearson val' + '\n')

            self.buildData_ic50()


    def buildData_ic50(self):
        """
        looks through file with outlier counts and saves it in a dictionary.
        """
        self.ccle = {}
        self.drug_list = {}

        with open(self.input1) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader)
            for row in tsvReader:
                a = row[0].split('-')
                b = a[3].replace('_', '-')
                self.ccle[b] = [row[1], row[2]]

        self.ccle_list =  self.ccle.keys()

        #print (self.ccle_list)

        with open(self.input2) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for row in tsvReader:
                #print(row)
                try:
                    self.drug_list[row[5]].append([row[3], row[9]])
                except KeyError:
                    self.drug_list[row[5]] = [[row[3], row[9]]]

        #print(self.drug_list)
        #print(len(self.drug_list['Erlotinib']))

        for key in self.drug_list:
            self.drug_analysis(key)



    def drug_analysis(self, key):
        """
        check overlap
        """
        tissue = {}

        for item in self.drug_list[key]:
            if item[0] in self.ccle_list:
                new_key_name =  key + '_' + self.ccle[item[0]][1]
                try:
                    tissue[new_key_name].append([item[1], self.ccle[item[0]][0], item[0]])
                except KeyError:
                    tissue[new_key_name] = [[item[1], self.ccle[item[0]][0], item[0]]]

        correlation_values = []

        # first get p-values for all drug-tissue combinations
        for key in tissue:
            if len(tissue[key]) >= 10:
                self.getSpearmanr(key, tissue[key], correlation_values)

        #print(correlation_values)
        # now get all corrected p-values

        self.multitestCorrection(correlation_values)



    def getSpearmanr(self, key, list, empty_list):

        x=[]
        y=[]

        for item in list:
            x.append(float(item[0]))
            y.append(float(item[1]))


        u, v = spearmanr(x,y)
        w, x = pearsonr(x,y)

        empty_list.append([key, u, v, w, x])





    def get_plots(self, xy_vals, corr_vals):

        plot_name = corr_vals[0]
        x=[]
        y=[]

        figure_width=4
        figure_height=4

        plt.figure(figsize=(figure_width,figure_height))

        panel_width=(1/figure_width) * 3
        panel_height=(1/figure_height) * 3



        panel1=plt.axes([0.15,0.1,panel_width,panel_height])

        sort_xy_vals = sorted(xy_vals, key = lambda x:float(x[0]))

        #print(sort_xy_vals)

        for item in sort_xy_vals:
            x.append(float(item[0]))
            y.append(float(item[1]))
            panel1.scatter(float(item[0]), float(item[1]), label=item[2],  marker = '.', s = 200, alpha = 0.6)

        panel1.text(.5,.7,'spearman: %.3f \n pval: %.3f \n pearson: %.3f \n pval: %.3f' % \
                    (corr_vals[1], corr_vals[2], corr_vals[3], corr_vals[4]), \
                    bbox={'facecolor':'w','pad':5}, \
                    ha="center", va="bottom", transform=plt.gca().transAxes)


        #select max and min values for tick_params

        x_max =  max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)


        panel1.set_xticks([x_min, x_max])
        panel1.set_yticks([y_min, y_max])


        panel1.plot(x, y, 'blue')
        panel1.set_title(plot_name)
        panel1.set_xlabel('log(IC50)')
        panel1.set_ylabel('outlier counts')
        #panel1.autoscale(enable=True, axis='y', tight=None)
        plt.savefig(plot_name, dpi=600)

        #plt.show()



    def multitestCorrection(self, list):
        final_values = []

        p_values1 = []
        p_values2 = []

        for item in list:
            p_values1.append(item[2])
            p_values2.append(item[4])

        a, b, c, d = multipletests(p_values1, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        e, f, g, h = multipletests(p_values2, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)



        for pvalue in range(0, len(b)):
            self.out_file.write(list[pvalue][0] +
                            '\t' + str(list[pvalue][1]) +
                            '\t' + str(list[pvalue][2]) +
                            '\t' + str(b[pvalue]) +
                            '\t' + str(list[pvalue][3]) +
                            '\t' + str(list[pvalue][4]) +
                            '\t' + str(f[pvalue]) +
                            '\n')



def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
