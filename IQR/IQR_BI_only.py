""" author: mcatamur@ucsc.edu """


"""
Sample Command: python IQR_tissue_outlier.py -i juncBase_table.txt/tsv
"""

#!/usr/bin/env python3


import csv, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
from collections import defaultdict



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
            self.parser.add_argument('-dPSI','--delta_PSI', type=int, action = 'store', help = 'delta PSI threshold', default=10.0)

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.dPSI = myCommandLine.args.delta_PSI
        self.findTissueIndex(self.inputFile)
        self.buildData_tissue(self.inputFile)
        self.plotEventsBar_allTissue()

    def findTissueIndex(self, inputFile):
        print("Building tissue index...")
        self.tissue_index = defaultdict(list)

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            self.header = next(tsvReader, None)
            for ccle in range(11,len(self.header)-1): #index of cell lines
                a = self.header[ccle].split("-") #split ccle tag
                if a[0] == 'BI':
                    # Hard code four letter codes
                    if a[1] == 'LUSC':
                        self.tissue_index['LUNG'].append(ccle)
                    if a[1] == 'LCLL':
                        self.tissue_index['LCLL'].append(ccle)
                    if a[1] == 'DLBC':
                        self.tissue_index['DLBC'].append(ccle)
                    if a[1] == 'MM':
                        self.tissue_index['MM'].append(ccle)
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


        #for item in self.tissue_index:
            #print (item, len(self.tissue_index[item]))

        self.total_tissue_outlier_counts = {key:0 for key, value in self.tissue_index.items()}

        #len(self.tissue_outlier_counts[key]) should equal len(self.tissue_index[key])

    def buildData_tissue(self, inputFile):
        print("Building IQR data..")
        #loop through the IQR method here
        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)
            row_num = 0
            for row in tsvReader:
                print(row_num)
                for tissue_type in self.tissue_index:
                    if len(self.tissue_index[tissue_type]) >= 7:
                        self.doIQR_part1(tissue_type, row)
                    else:
                        pass
                row_num += 1


    def doIQR_part1(self, tissue_type, row):
        psi_values = []
        for index in self.tissue_index[tissue_type]:
            try:
                psi_values.append(float(row[index]))
            except ValueError:
                continue


        psi_values.sort(key = int)
        row_len = len(psi_values)

        if (row_len % 2) != 0 and row_len > 7:
            q2 = psi_values[(row_len//2) + 1]
            q1 = (psi_values[(row_len//2)//2] + psi_values[((row_len//2)//2)-1]) / 2
            q3 = (psi_values[(row_len//2) + ((row_len//2)//2)] + psi_values[(row_len//2) + (((row_len//2)//2)+1)]) / 2
            self.doIQR_part2(psi_values, tissue_type, q1, q2, q3, row)

        elif row_len == 7: #Special case
            q2 = psi_values[(row_len//2) + 1]
            q1 = psi_values[(row_len//2)//2]
            q3 = psi_values[(row_len//2)+(((row_len//2)//2)+1)]
            self.doIQR_part2(psi_values, tissue_type, q1, q2, q3, row)

        elif (row_len % 2) == 0 and row_len >= 8:
            q2 = (psi_values[row_len//2] + psi_values[(row_len//2) - 1]) / 2
            if (row_len//2) % 2 == 0:
                q1 = psi_values[((row_len//2)//2) - 1]
                q3 = psi_values[(row_len//2) + (((row_len//2)//2)+1)]
                self.doIQR_part2(psi_values, tissue_type, q1, q2, q3, row)
            elif (row_len//2) % 2 != 0:
                q1 = psi_values[(row_len//2)//2]
                q3 = psi_values[(row_len//2)+(((row_len//2)//2)+1)]
                self.doIQR_part2(psi_values, tissue_type, q1, q2, q3, row)

        elif row_len <= 6:
            pass


    def doIQR_part2(self, psi_values, tissue, q1, q2, q3, row):

        iqr = q3-q1

        #instead of indexing through the psi_values, let's go by the tissue index instead
        #that way we can also acces to which cell line the value belonged to

        index =  self.tissue_index[tissue]

        for i in range(11,len(self.header)-1):
            if i in index:
                try:
                    if float(row[i]) < (q1-1.5*iqr) or float(row[i]) > (q3+1.5*iqr):
                        if (q2-float(row[i])) > self.dPSI or (float(row[i])-q2) > self.dPSI:
                            self.total_tissue_outlier_counts[tissue] += 1
                        else:
                            pass
                    else:
                        pass
                except ValueError:
                    pass






    def plotEventsBar_allTissue(self):
        print(self.total_tissue_outlier_counts)
        print(self.tissue_index['BREAST'])
        self.colors = ['#FFD54F', '#FFF176', '#81C784', '#4DD0E1', '#64B5F6', '#9575CD', '#F06292', '#EF5350', '#FF8A65', '#A1887F','#90A4AE']
        tuple_tissue_events = list(self.total_tissue_outlier_counts.items())
        for tissue in tuple_tissue_events: #(event, counts)
            plt.bar((tuple_tissue_events.index(tissue)), height = tissue[1], color = self.colors[6], label = tissue[0] )
            #plt.bar((x+1.5), height = f_events[x][1], color = self.colors[x], label = f_events[x][0])
            #print (event)
            #print (f_events.index(event))
            #print (self.filteredEvents[event[1]])


        plt.xticks(range(len(self.total_tissue_outlier_counts)), list(self.total_tissue_outlier_counts.keys()), rotation='vertical')
        plt.title('Outlier Splicing Events across all tissue types with dPSI > %s' % (str(self.dPSI)))
        plt.ylabel('Outlier Splicing Counts')
        figname_events_bar = 'all_tissue_types_events_bar_IQR_adjust.png'
        plt.savefig(figname_events_bar)
        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
