""" author: mcatamur@ucsc.edu """


"""
Sample Command: python splicing_mod_tissue.py -i juncBase_table.txt/tsv
"""

#!/usr/bin/env python3


import csv, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import seaborn as sns; sns.set()
from sklearn.decomposition import PCA
from matplotlib.font_manager import FontProperties


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

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.findTissueIndex(self.inputFile)
        self.buildData_tissue(self.inputFile)
        self.adjust_NA_values()
        self.remove_NA_values()
        self.Zscore()
        self.doPCA()


    def findTissueIndex(self, inputFile):

        self.tissue_index = defaultdict(list)

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for ccle in range(11,len(header)-1): #index of cell lines
                a = header[ccle].split("-") #split ccle tag
                if a[0] != 'BI':
                    d = a[0].split("_")
                    del d[0] #remove fh
                    del d[0] #remove 2nd code
                    tissue_name = '_'.join(d)
                    self.tissue_index[tissue_name].append(ccle)
                else:
                    pass



    def buildData_tissue(self, inputFile):
        self.raw_PCA_data = defaultdict(list)

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)

            array_length = (len(header)-1)-11

            row_num = 1
            for row in tsvReader:
                for ccle in range(11,len(header)-1):
                    try:
                        self.raw_PCA_data[row_num].append(float(row[ccle]))
                    except ValueError:
                        self.raw_PCA_data[row_num].append(row[ccle])
                row_num += 1



    def adjust_NA_values(self):

        for event in self.raw_PCA_data:
            NA_count = 0
            for value in self.raw_PCA_data[event]:
                if value == 'NA':
                    NA_count += 1

            #Check if NaNs are only 1% of the data.
            if float(NA_count/len(self.raw_PCA_data[event])) <= 0.01:
                #If 1%, change all NA values to the median
                #Append all numerical data to a list, then sort
                psi_values = []
                for value in self.raw_PCA_data[event]:
                    if value != 'NA':
                        psi_values.append(value)

                sorted_psi_values = sorted(psi_values)

                #determine the median
                row_len = len(sorted_psi_values)
                if row_len % 2 == 0: #even set
                    median = (sorted_psi_values[row_len//2] + sorted_psi_values[(row_len//2)-1]) / 2

                if row_len % 2 != 0: #odd set
                    median = sorted_psi_values[row_len//2]

                #change all NAs to the median value

                for i in range(0, row_len-1):
                    if self.raw_PCA_data[event][i] == 'NA':
                        self.raw_PCA_data[event][i] = median
            else:
                pass


    def remove_NA_values(self):
            self.adjusted_PCA_data = {}

            for event in self.raw_PCA_data:
                if 'NA' not in self.raw_PCA_data[event]:
                    self.adjusted_PCA_data[event] = self.raw_PCA_data[event]



    def Zscore(self):
        for event in self.adjusted_PCA_data:
            mean = sum(self.adjusted_PCA_data[event])/len(self.adjusted_PCA_data[event])

            #calculate variance
            variance = 0
            for value in self.adjusted_PCA_data[event]:
                variance += (value - mean)**2
                variance = variance / (len(self.adjusted_PCA_data[event])-1)

            #find standard deviation
            std_dev = variance**(1/2)

            #Now, transform PSI values into Z scores

            for i in range(0, len(self.adjusted_PCA_data[event])-1):
                self.adjusted_PCA_data[event][i] = ((self.adjusted_PCA_data[event][i]) - mean) / std_dev



        #Last step is to modify the dictionary from events to individual CCLEs
        self.PCA_data_tissue = {}


        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            #index for PSI values
            i = 0
            for ccle in range(11,len(header)-1):
                if i < len(self.raw_PCA_data.keys()):
                    for event in self.adjusted_PCA_data:
                        try:
                            self.PCA_data_tissue[header[ccle]].append(self.adjusted_PCA_data[event][i])
                        except KeyError:
                            self.PCA_data_tissue[header[ccle]] = [self.adjusted_PCA_data[event][i]]
                    i += 1



    def doPCA(self):
        colors = {'CENTRAL_NERVOUS_SYSTEM': '#34495E', \
                   'AUTONOMIC_GANGLIA' : '#808B96', \
                   'BONE_AND_SOFT_TISSUE' : '#1A5276', \
                   'SOFT_TISSUE' : '#3498DB', \
                   'BONE' : '#85C1E9', \
                   'CERVICAL' : '#2ECC71', \
                   'ENDOMETRIUM' : '#27AE60', \
                   'OVARY' : '#1ABC9C', \
                   'PROSTATE' : '#16A085', \
                   'URINARY_TRACT' : '#D35400', \
                   'KIDNEY' : '#E67E22', \
                   'BILIARY' : '#F39C12', \
                   'PANCREAS': '#F1C40F', \
                   'THYROID' :'#873600', \
                   'LIVER' : '#E74C3C', \
                   'LARGE_INTESTINE' : '#A93226', \
                   'BREAST' : '#6C3483', \
                   'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE' : '#A569BD', \
                   'OESOPHAGUS' : '#784212', \
                   'UPPER_AERODIGESTIVE_TRACT' : '#9C640C', \
                   'LUNG' : '#EC7063', \
                   'PLEURA' : '#F1948A', \
                   'MESOTHELIOMA' : '#17202A', \
                   'SKIN' : '#626567', \
                   'STOMACH' : '#7B241C'}

        #pca = PCA(n_components=2)
        X_raw = np.asarray((list(self.PCA_data_tissue.values())))
        X = PCA(n_components=2).fit_transform(X_raw)
        #pca.fit(X)
        #X = pca.fit_transform(X)
        #X = pca.fit_transform(X)
        print (X.shape)
        #print (X)



        #print one of each for legend
        for tissue in self.tissue_index:
            i = 0
            if i >= 0:
                plt.scatter(X[self.tissue_index[tissue][i]][0], X[self.tissue_index[tissue][i]][1], \
                            label = tissue, color = colors[tissue], marker = '.', s = 60, alpha = 0.5)



        for tissue in self.tissue_index:
            for index in self.tissue_index[tissue]:
                plt.scatter(X[index-11][0], X[index-11][1], color = colors[tissue], marker = '.', s = 60, alpha = 0.5)


        legend = list(colors.keys())
        plt.autoscale(enable=True, axis='both', tight=None)
        #plt.axis([-20, 80, -40, 85])
        plt.ylabel('PC 1')
        plt.xlabel('PC 2')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()



def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
