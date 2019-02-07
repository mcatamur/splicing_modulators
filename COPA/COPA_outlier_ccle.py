""" author: mcatamur@ucsc.edu """


"""
Sample Command: python splicing_mod.py -i CCLE..
"""
#!/usr/bin/env python3


import csv, time, argparse
start_time = time.time()
import matplotlib
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import numpy as np
import math


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
            self.parser.add_argument('-c','--ccle', type=str, action = 'store', help = 'cell line to analyze', default='BI-LGG-CCLE-MHH_NB_11-Tumor')
            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class outputAnalyzer() :

    def __init__(self, myCommandLine):
        self.inputFile = myCommandLine.args.inputFile
        self.ccle = myCommandLine.args.ccle
        self.buildData(self.inputFile)
        self.plotEventsBar()
        self.plotEventsBarCCLE()
        #self.plotScatterAll()
        #self.plotScatterEvent()

    def buildData(self, inputFile):
        """
        Build dictionaries in this method for stacked bar plot and volcano scatter plot
        """
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


    def plotEventsBar(self):
        """
        This method creates a bar graph comparison between all significant alternative splicing events
        across all cell lines and with PSI > 0 for the --ccle input
        """
        #Display number of events as a fraction of a whole
        totalEvents = sum(self.allEvents.values())
        totalAdjusted = sum(self.filteredEvents.values())

        #print (totalAdjusted)
        #print (totalEvents)

        frac_dict_all = {key:0 for key,value in self.allEvents.items()}
        frac_dict_filtered = {key:0 for key,value in self.allEvents.items()}

        for event in frac_dict_all:
            frac_dict_all[event] = (self.allEvents[event]/totalEvents) * 100
            frac_dict_filtered[event] = (self.filteredEvents[event]/totalAdjusted) * 100



        #aesthetics
        xlabel = ['All AS event counts \n across all \n cell lines',
                   'AS event counts for \n %s \n with PSI > 0' % (self.ccle)]

        self.colors = ['#FFD54F', '#FFF176', '#81C784', '#4DD0E1', '#64B5F6', '#9575CD', '#F06292', '#EF5350', '#FF8A65', '#A1887F','#90A4AE']

        #plot the results in a stacked bar graph


        #Turn the dict to a tuple. That way it is ordered and is subscriptable.
        a_events = list(frac_dict_all.items())
        self.f_events = list(frac_dict_filtered.items())
        x_vals = [0.5, 2.5] #distance of bars from each other

        #Plot the Top bar first
        plt.bar(x_vals[0], a_events[0][1], bottom = 0, color = self.colors[0], label = a_events[0][0])
        plt.bar(x_vals[1], self.f_events[0][1], bottom = 0, color = self.colors[0]) #same label for both
        a_height = [a_events[0][1]]
        f_height = [self.f_events[0][1]]

        #Plot the rest
        for x in range(1, len(frac_dict_all.keys())):
            try:
                plt.bar(x_vals[0], height = a_events[x][1], bottom=sum(a_height), color=self.colors[x], label=a_events[x][0])
                a_height.append(a_events[x][1])
                plt.bar(x_vals[1], height = self.f_events[x][1], bottom=sum(f_height), color=self.colors[x])
                f_height.append(self.f_events[x][1])
            except IndexError:
                continue

        figname_bar = self.ccle + '_' + 'bar' + '.png'
        plt.xticks(x_vals, xlabel)
        plt.ylabel('Proportion')
        plt.legend()
        plt.axis([0, 5.5, 0, 100])
        plt.title('Signinficant Alternative Splicing Event Counts')
        plt.savefig(figname_bar)
        plt.show()



    def plotEventsBarCCLE(self):
        for event in self.f_events: #(event, counts)
            plt.bar((self.f_events.index(event)), height = self.filteredEvents[event[0]], color = self.colors[self.f_events.index(event)], label = event[0] )
            #plt.bar((x+1.5), height = self.f_events[x][1], color = self.colors[x], label = self.f_events[x][0])
            #print (event)
            #print (self.f_events.index(event))
            #print (self.filteredEvents[event[1]])


        plt.xticks(range(len(self.filteredEvents)), list(self.filteredEvents.keys()), rotation='vertical')
        plt.title('Signinficant Alternative Splicing Event Counts \n for %s with PSI > 0' % (self.ccle))
        plt.ylabel('Signinficant AS Event Counts')
        figname_events_bar = self.ccle + '_' + 'events' + '_' + 'bar' + '.png'
        plt.savefig(figname_events_bar)
        plt.show()





    def plotOneEventScatter(self, event):
        y_all = []
        x_all = []

        for valueSet in self.oneEvent_all[event]:
            x_all.append(valueSet[0])
            y_all.append(-(math.log(valueSet[1])))

        #print (self.oneEvent_filtered)

        y_filtered = []
        x_filtered = []

        for valueSet in self.oneEvent_filtered[event]:
            x_filtered.append(valueSet[0])
            y_filtered.append(-(math.log(valueSet[1])))

        legend = ['All other \'%s\' events.\n %s total %s events' % (event, len(x_all), event),
                  '\'%s\' events with \n corrected p-value above \n the threshold. p < %s. \n %s total events' % (event, self.pVal, len(x_filtered))]

        plt.scatter(x_filtered,y_filtered, label = legend[1], color = self.colors[5], marker = '.', s = 60, alpha = 0.6)
        plt.scatter(x_all,y_all, label = legend[0], color = self.colors[1], marker = '.', s = 60, alpha = 0.5)

        plt.autoscale(enable=True, axis='both', tight=None)
        plt.ylabel('-log(corrected p-value)')
        plt.xlabel('Delta PSI Value')
        plt.legend()



    def plotAllEventsScatter(self):

        y_all = []
        x_all = []

        for value in self.all_pVal:
            x_all.append(value[0])
            y_all.append(-(math.log(value[1]))) #Modify the corrected p-values to -log(corrected p-value)

        y_filtered = []
        x_filtered = []

        #print (self.filtered_pVal)

        for value in self.filtered_pVal:
            x_filtered.append(value[0])
            y_filtered.append(-(math.log(value[1])))


        plt.scatter(x_filtered,y_filtered, label = self.legend[1], color = self.colors[5], marker = '.', s = 60, alpha = 0.6)
        plt.scatter(x_all,y_all, label = self.legend[0], color = self.colors[1], marker = '.', s = 60, alpha = 0.5)

        plt.autoscale(enable=True, axis='both', tight=None)
        #plt.axis([-50, 50, 0, 20])
        plt.ylabel('-log(corrected p-value)')
        plt.xlabel('Delta PSI Value')
        plt.legend()

        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myCommandLine = CommandLine(myCommandLine)
    myFileReader = outputAnalyzer(myCommandLine)
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
