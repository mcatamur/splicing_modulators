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
        self.buildData_tissue(self.inputFile)
        self.plotEventsBar()
        self.plotEventsBar_allTissue()
        self.plotEventsBar_oneTissue()


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

        #for item in self.tissue_index:
            #print (item, len(self.tissue_index[item]))

        self.tissue_event_counts = {key:0 for key, value in self.tissue_index.items()}


    def buildData_tissue(self, inputFile):
        self.allEvents_all_tissue_outlier = {}

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)

        #first we have to create the keys in the dictionary.
        #take into account if count is 0
            for row in tsvReader:
                self.allEvents_all_tissue_outlier[row[1]] = 0 #all events, regardless of tissue type or outlier values

        self.allEvents_spec_tissue_outlier = {key:0 for key,value in self.allEvents_all_tissue_outlier.items()}
        self.events_spec_tissue = {key:0 for key,value in self.allEvents_all_tissue_outlier.items()}

        with open(self.inputFile) as file:
            tsvReader = csv.reader(file, dialect = 'excel-tab')
            next(tsvReader, None)
            for row in tsvReader:
                self.doIQR_part1(row)


    def doIQR_part1(self, row):
        psi_values = []
        for tissue in self.tissue_index:
            if len(self.tissue_index[tissue]) >= 7:
                for index in self.tissue_index[tissue]:
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
            self.doIQR_part2(q1, q2, q3, row)

        elif row_len == 7: #Special case
            q2 = psi_values[(row_len//2) + 1]
            q1 = psi_values[(row_len//2)//2]
            q3 = psi_values[(row_len//2)+(((row_len//2)//2)+1)]
            self.doIQR_part2(q1, q2, q3, row)

        elif (row_len % 2) == 0 and row_len >= 8:
            q2 = (psi_values[row_len//2] + psi_values[(row_len//2) - 1]) / 2
            if (row_len//2) % 2 == 0:
                q1 = psi_values[((row_len//2)//2) - 1]
                q3 = psi_values[(row_len//2) + (((row_len//2)//2)+1)]
                self.doIQR_part2(q1, q2, q3, row)
            elif (row_len//2) % 2 != 0:
                q1 = psi_values[(row_len//2)//2]
                q3 = psi_values[(row_len//2)+(((row_len//2)//2)+1)]
                self.doIQR_part2(q1, q2, q3, row)

        elif row_len <= 6:
            pass

        #if row_len > 6:
            #print (q1)
            #print(q3)

    def doIQR_part2(self, q1, q2, q3, row):

        iqr = q1-q3

        for key in self.tissue_index:
            if len(self.tissue_index[key]) >= 7:
                for index in self.tissue_index[key]:
                    try:
                        if float(row[index]) < (q1-1.5*iqr) or float(row[index]) > (q3+1.5*iqr):
                            if (q2 - float(row[index])) > self.dPSI or (float(row[index])-q2) > self.dPSI:
                                self.allEvents_all_tissue_outlier[row[1]] += 1
                                self.tissue_event_counts[key] += 1
                    except ValueError:
                        continue

            if key == self.tissue:
                for index in self.tissue_index[key]:
                    try:
                        if float(row[index]) < (q1-1.5*iqr) or float(row[index]) > (q3-1.5*iqr):
                            if (q1 - float(row[index])) > self.dPSI or (float(row[index])-q3) > self.dPSI:
                                self.allEvents_spec_tissue_outlier[row[1]] += 1
                                self.events_spec_tissue[row[1]] += 1
                    except ValueError:
                        continue



    def plotEventsBar(self):
        """
        This method creates a bar graph comparison between all significant alternative splicing events
        across all cell lines and with dPSI > self.dPSI in IQR Analysis
        """


        total_all_tissues =  sum(self.allEvents_all_tissue_outlier.values())
        total_spec_tissue = sum(self.allEvents_spec_tissue_outlier.values())

        frac_spec_tissue = {key:0 for key,value in self.allEvents_all_tissue_outlier.items()}
        frac_all_tissue = {key:0 for key,value in self.allEvents_all_tissue_outlier.items()}

        event_keys = self.allEvents_all_tissue_outlier.keys()
        for event in event_keys:
            frac_all_tissue[event] = (float(self.allEvents_all_tissue_outlier[event]/total_all_tissues)) * 100.0
            frac_spec_tissue[event] = (float(self.allEvents_spec_tissue_outlier[event]/total_spec_tissue)) * 100.0


        #aesthetics
        xlabel = ['Outlier splicing counts \n across all tissue types \n with dPSI > %s' % (self.dPSI),
                   'Outlier splicing counts for \n %s \n with dPSI > %s' % (self.tissue, self.dPSI)]

        self.colors = ['#FFD54F', '#FFF176', '#81C784', '#4DD0E1', '#64B5F6', '#9575CD', '#F06292', '#EF5350', '#FF8A65', '#A1887F','#90A4AE']

        #plot the results in a stacked bar graph


        #Turn the dict to a tuple. That way it is ordered and is subscriptable.
        a_events = list(frac_all_tissue.items())
        f_events = list(frac_spec_tissue.items())
        x_vals = [0.5, 2.5] #distance of bars from each other



        #Plot the Top bar first
        plt.bar(x_vals[0], a_events[0][1], bottom = 0, color = self.colors[0], label = a_events[0][0])
        plt.bar(x_vals[1], f_events[0][1], bottom = 0, color = self.colors[0]) #same label for both
        a_height = [a_events[0][1]]
        f_height = [f_events[0][1]]

        #Plot the rest
        for x in range(1, len(frac_all_tissue.keys())):
            try:
                plt.bar(x_vals[0],height = a_events[x][1], bottom=sum(a_height), color=self.colors[x], label=a_events[x][0])
                a_height.append(a_events[x][1])
                plt.bar(x_vals[1], height = f_events[x][1], bottom=sum(f_height), color=self.colors[x])
                f_height.append(f_events[x][1])
            except IndexError:
                continue

        figname_bar = self.tissue + '_' + 'bar' + '.png'
        plt.xticks(x_vals, xlabel)
        plt.ylabel('Proportion')
        plt.legend()
        plt.axis([0, 5.5, 0, 100])
        plt.title('Outlier Splicing Counts')
        plt.savefig(figname_bar)
        plt.show()



    def plotEventsBar_allTissue(self):
        tuple_tissue_events = list(self.tissue_event_counts.items())
        for tissue in tuple_tissue_events: #(event, counts)
            plt.bar((tuple_tissue_events.index(tissue)), height = tissue[1], color = self.colors[4], label = tissue[0] )
            #plt.bar((x+1.5), height = f_events[x][1], color = self.colors[x], label = f_events[x][0])
            #print (event)
            #print (f_events.index(event))
            #print (self.filteredEvents[event[1]])


        plt.xticks(range(len(self.tissue_event_counts)), list(self.tissue_event_counts.keys()), rotation='vertical')
        plt.title('Outlier Splicing Events across all tissue types with dPSI > %s' % (str(self.dPSI)))
        plt.ylabel('Outlier Splicing Counts')
        figname_events_bar = 'all_tissue_types_events_bar.png'
        plt.savefig(figname_events_bar)
        plt.show()



    def plotEventsBar_oneTissue(self):
        tuple_spec_tissue_events = list(self.events_spec_tissue.items())
        for tissue in tuple_spec_tissue_events: #(event, counts)
            plt.bar((tuple_spec_tissue_events.index(tissue)), height = tissue[1], color = self.colors[tuple_spec_tissue_events.index(tissue)], label = tissue[0] )

        plt.xticks(range(len(self.events_spec_tissue)), list(self.events_spec_tissue.keys()), rotation='vertical')
        plt.title('Outlier Splicing Events for %s with PSI > %s' % (self.tissue, str(self.dPSI)))
        plt.ylabel('Outlier Splicing Event Counts')
        figname_events_bar = self.tissue + '_' + 'events' + '_' + 'bar' + '.png'
        plt.savefig(figname_events_bar)
        plt.show()
        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
