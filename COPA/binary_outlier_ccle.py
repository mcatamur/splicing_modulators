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
            self.parser.add_argument('-o1', '--outputFile1', action = 'store', help = 'first output file')
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
        self.outputFile1 = myCommandLine.args.outputFile1
        self.CCLE_dict = {}
        self.CCLE_gene_dict = {}
        self.buildData(self.inputFile)
        self.plotEventsBar_allTissue()





    def buildData(self, inputFile):
        """
        Build dictionaries in this method for stacked bar plot and volcano scatter plot
        """

        self.tissue_index = defaultdict(list)

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            #first we have to create the keys in the dictionary or each cancer cell line
            header = next(tsvReader)
            for ccle in header:
                if ccle != 'gene':
                    a = ccle.split("-") #split ccle tag
                    # Hard code four letter codes
                    if a[0] == 'BI':
                        if a[1] == 'LUSC':
                            self.tissue_index['LUNG'].append(header.index(ccle))
                        if a[1] == 'LCLL' or a[1] == 'DLBC' or a[1] == 'MM':
                            self.tissue_index['HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'].append(header.index(ccle))
                        if a[1] == 'LGG':
                            self.tissue_index['CENTRAL_NERVOUS_SYSTEM'].append(header.index(ccle))
                        if a[1] == 'BRCA':
                            self.tissue_index['BREAST'].append(header.index(ccle))
                        if a[1] == 'COAD':
                            self.tissue_index['LARGE_INTESTINE'].append(header.index(ccle))
                        if a[1] == 'SARC':
                            self.tissue_index['BONE_AND_SOFT_TISSUE'].append(header.index(ccle))
                        if a[1] == 'OV':
                            self.tissue_index['OVARY'].append(header.index(ccle))
                        if a[1] == 'SKCM':
                            self.tissue_index['SKIN'].append(header.index(ccle))
                        if a[1] == 'PAAD':
                            self.tissue_index['PANCREAS'].append(header.index(ccle))
                        if a[1] == 'HNSC':
                            self.tissue_index['UPPER_AERODIGESTIVE_TRACT'].append(header.index(ccle))
                        if a[1] == 'KIRC':
                            self.tissue_index['KIDNEY'].append(header.index(ccle))
                        if a[1] == 'CESC':
                            self.tissue_index['CERVICAL'].append(header.index(ccle))
                        if a[1] == 'BLCA':
                            self.tissue_index['URINARY_TRACT'].append(header.index(ccle))
                        if a[1] == 'ESCA':
                            self.tissue_index['OESOPHAGUS'].append(header.index(ccle))
                        if a[1] == 'LIHC':
                            self.tissue_index['BILIARY'].append(header.index(ccle))
                        if a[1] == 'STAD':
                            self.tissue_index['STOMACH'].append(header.index(ccle))
                        if a[1] == 'THCA':
                            self.tissue_index['THYROID'].append(header.index(ccle))
                        if a[1] == 'PRAD':
                            self.tissue_index['PROSTATE'].append(header.index(ccle))
                        if a[1] == 'MESO':
                            self.tissue_index['MESOTHELIOMA'].append(header.index(ccle))

                    else:
                        d = a[0].split("_")
                        del d[0] #remove fh
                        del d[0] #remove 2nd code
                        tissue_name = '_'.join(d)
                        self.tissue_index[tissue_name].append(header.index(ccle))


        #Now we look through the binary table, if 1: +=1,
        total_rows = 0
        self.tissue_outlier_counts = {key:0 for key,value in self.tissue_index.items()}

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            header = next(tsvReader)
            for row in tsvReader:
                for tissue in self.tissue_index:
                    for index in self.tissue_index[tissue]:
                        if int(row[index]) == 1:
                            self.tissue_outlier_counts[tissue] += 1

        #sorted keys
        mean_ccle_dict = sum(self.tissue_outlier_counts.values())/len(self.tissue_outlier_counts)


        with open(self.outputFile1, mode ='w') as out_file1:
            tsvWriter = csv.writer(out_file1, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            tsvWriter.writerow(['mean:', mean_ccle_dict])
            sorted_dict = sorted(self.tissue_outlier_counts.items(), key=operator.itemgetter(1))
            for key in sorted_dict:
                tsvWriter.writerow([key[0], key[1]])


    def plotEventsBar_allTissue(self):
        colors = ['#FFD54F', '#FFF176', '#81C784', '#4DD0E1', '#64B5F6', '#9575CD', '#F06292', '#EF5350', '#FF8A65', '#A1887F','#90A4AE']

        tuple_tissue_events = list(self.tissue_outlier_counts.items())
        for tissue in tuple_tissue_events: #(event, counts)
            plt.bar((tuple_tissue_events.index(tissue)), height = tissue[1], color = colors[2], label = tissue[0] )

        plt.xticks(range(len(self.tissue_outlier_counts)), list(self.tissue_outlier_counts.keys()), rotation='vertical')
        plt.title('Outlier Splicing Events across all tissue types')
        plt.ylabel('Outlier Splicing Counts')
        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = binaryOutliers(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
