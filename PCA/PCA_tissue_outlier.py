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
        self.doPCA()


    def findTissueIndex(self, inputFile):

        self.tissue_index = defaultdict(list)

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            for ccle in range(11,len(header)-1): #index of cell lines
                a = header[ccle].split("-") #split ccle tag

                # Hard code four letter codes
                if a[0] == 'BI':
                    if a[1] == 'LUSC':
                        self.tissue_index['LUNG'].append(ccle)
                    if a[1] == 'LCLL' or a[1] == 'DLBC' or a[1] == 'MM':
                        self.tissue_index['HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'].append(ccle)
                    if a[1] == 'LGG':
                        self.tissue_index['CENTRAL_NERVOUS_SYSTEM'].append(ccle)
                    if a[1] == 'BRCA':
                        self.tissue_index['BREAST'].append(ccle)
                    if a[1] == 'COAD':
                        self.tissue_index['LARGE_INTESTINE'].append(ccle)
                    if a[1] == 'SARC':
                        self.tissue_index['BONE_AND_SOFT_TISSUE'].append(ccle)
                    if a[1] == 'OV':
                        self.tissue_index['OVARY'].append(ccle)
                    if a[1] == 'SKCM':
                        self.tissue_index['SKIN'].append(ccle)
                    if a[1] == 'PAAD':
                        self.tissue_index['PANCREAS'].append(ccle)
                    if a[1] == 'HNSC':
                        self.tissue_index['UPPER_AERODIGESTIVE_TRACT'].append(ccle)
                    if a[1] == 'KIRC':
                        self.tissue_index['KIDNEY'].append(ccle)
                    if a[1] == 'CESC':
                        self.tissue_index['CERVICAL'].append(ccle)
                    if a[1] == 'BLCA':
                        self.tissue_index['URINARY_TRACT'].append(ccle)
                    if a[1] == 'ESCA':
                        self.tissue_index['OESOPHAGUS'].append(ccle)
                    if a[1] == 'LIHC':
                        self.tissue_index['BILIARY'].append(ccle)
                    if a[1] == 'STAD':
                        self.tissue_index['STOMACH'].append(ccle)
                    if a[1] == 'THCA':
                        self.tissue_index['THYROID'].append(ccle)
                    if a[1] == 'PRAD':
                        self.tissue_index['PROSTATE'].append(ccle)
                    if a[1] == 'MESO':
                        self.tissue_index['MESOTHELIOMA'].append(ccle)

                else:
                    pass



    def buildData_tissue(self, inputFile):
        self.raw_PCA_data = defaultdict(list)
        self.reduced_PCA_data = {}

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
                        pass
                row_num += 1

        self.checkArraySize(array_length)


    def checkArraySize(self, array_length):
        for key in self.raw_PCA_data:
            if len(self.raw_PCA_data[key]) == array_length:
                self.reduced_PCA_data[key] = self.raw_PCA_data[key]

        #modify the dictionary from events to tissue type
        self.PCA_data_tissue = {}


        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)
            #index for PSI values
            i = 0
            for ccle in range(11,len(header)-1):
                if i < len(self.raw_PCA_data.keys()):
                    for event in self.reduced_PCA_data:
                        try:
                            self.PCA_data_tissue[header[ccle]].append(self.reduced_PCA_data[event][i])
                        except KeyError:
                            self.PCA_data_tissue[header[ccle]] = [self.reduced_PCA_data[event][i]]
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

        pca = PCA(n_components=2)
        X = np.asarray((list(self.PCA_data_tissue.values())))
        X = pca.fit_transform(X)
        print (X.shape)




        #print one of each for legend
        for tissue in self.tissue_index:
            i = 0
            if tissue == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE' and i >= 0:
                plt.scatter(X[self.tissue_index[tissue][i]][0], X[self.tissue_index[tissue][i]][1], \
                            label = tissue, color = colors[tissue], marker = '.', s = 60, alpha = 0.5)



        for tissue in self.tissue_index:
            if tissue == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE':
                for index in self.tissue_index[tissue]:
                    plt.scatter(X[index-11][0], X[index-11][1], color = colors[tissue], marker = '.', s = 60, alpha = 0.5)


        legend = list(colors.keys())
        #plt.autoscale(enable=True, axis='both', tight=None)
        plt.axis([-20, 80, -40, 85])
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
