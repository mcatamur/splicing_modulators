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
from sklearn.manifold import TSNE
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
            self.parser.add_argument('-o', '--outputFile', action = 'store', help = 'outputfile')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.outputFile = myCommandLine.args.outputFile
        self.findTissueIndex(self.inputFile)
        #self.buildData_tissue(self.inputFile)
        self.createTuples(self.inputFile)
        #self.adjust_NA_values()
        self.remove_NA_values()
        #self.Zscore()
        self.modifyFinalArray()
        self.dotSNE()


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
                    d = a[0].split("_")
                    del d[0] #remove fh
                    del d[0] #remove 2nd code
                    tissue_name = '_'.join(d)
                    self.tissue_index[tissue_name].append(ccle)

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

    def createTuples(self, inputFile):
        #modify each row first

        #create a list of lists for the events
        self.list_of_all_events = []

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader, None)

            for row in tsvReader: #for every row
                current_row = []
                for ccle in range(11,len(header)-1):
                    try:
                        current_row.append(float(row[ccle]))
                    except ValueError:
                        current_row.append(row[ccle])

                self.adjust_NA_values(current_row)
        print("Finished adjusting non-PSI values...")




    def adjust_NA_values(self, current_row):

        NA_count = 0
        for value in current_row:
            if value == 'NA':
                NA_count += 1

        #Check if NaNs are only 1% of the data.
        if float(NA_count/len(current_row)) <= 0.01:
            #If 1%, change all NA values to the median
            #Append all numerical data to a list, then sort
            psi_values = []
            for value in current_row:
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
                if current_row[i] == 'NA':
                    current_row[i] = median


        else:
            pass


        #now append this row to self.list_of_all_events
        self.list_of_all_events.append(current_row)




    def remove_NA_values(self):
        self.list_of_all_events_2 = []
        print("Removing non-PSI values...")
        for row in self.list_of_all_events:
            if 'NA' not in row:
                self.list_of_all_events_2.append(row)



    def Zscore(self):

        print("Converting to Z-score values...")
        for row in self.list_of_all_events_2:
            mean = sum(row)/len(row)

            #calculate variance
            variance = 0
            for value in row:
                variance += (value - mean)**2
                variance = variance / (len(row)-1)

            #find standard deviation
            std_dev = variance**(1/2)

            #Now, transform PSI values into Z scores

            for i in range(0, len(row)-1):
                row[i] = (row[i] - mean) / std_dev


    def modifyFinalArray(self):
        print("Modifying final array...")
        #Last step is to modify the dictionary from events to individual CCLEs
        self.final_tSNE_data = []

        row_len = len(self.list_of_all_events_2[0])
        for i in range(0, row_len-1):
            i_list = []
            for row in self.list_of_all_events_2:
                i_list.append(row[i])
            self.final_tSNE_data.append(i_list)

        #check point stuff:
        print(len(self.final_tSNE_data))
        #should be the same number
        print(len(self.final_tSNE_data[1]))
        print(len(self.final_tSNE_data[444]))
        print(len(self.final_tSNE_data[666]))
        print(len(self.final_tSNE_data[738]))



    def dotSNE(self):

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
        X_raw = np.asarray(self.final_tSNE_data)
        print(self.final_tSNE_data[444][1])
        print(X_raw[444][1])
        X = TSNE(n_components=2).fit_transform(X_raw)
        print(X)
        #pca.fit(X)
        #X = pca.fit_transform(X)
        #X = pca.fit_transform(X)
        #print (X.shape)
        #print (X)

        with open(self.outputFile, mode ='w') as out_file1:
            tsvWriter = csv.writer(out_file1, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for x in X:
                tsvWriter.writerow(["(", [x[0], ",", x[1]], ")"])



        #print one of each for legend
        for tissue in self.tissue_index:
            i = 0
            if i >= 0:
                plt.scatter(X[self.tissue_index[tissue][i]][0], X[self.tissue_index[tissue][i]][1], \
                            label = tissue, color = colors[tissue], marker = '.', s = 60, alpha = 0.5)



        for tissue in self.tissue_index:
            for index in self.tissue_index[tissue]:
                try:
                    plt.scatter(X[index-11][0], X[index-11][1], color = colors[tissue], marker = '.', s = 60, alpha = 0.5)
                except IndexError:
                    pass

        legend = list(colors.keys())
        plt.autoscale(enable=True, axis='both', tight=None)
        #plt.axis([-20, 80, -40, 85])
        plt.ylabel('x[1]')
        plt.xlabel('x[2]')
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
