""" author: mcatamur@ucsc.edu """


"""
Sample Command: python IQR_tissue_outlier.py -i juncBase_table.txt/tsv
"""

#!/usr/bin/env python3


import csv, time, argparse
start_time = time.time()
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np
import pandas as pd


class outputAnalyzer() :

    def __init__(self):
        self.plotSwarm()


    def plotSwarm(self):
        ccle_df= pd.read_csv('swarm_plot_data_positive_control.txt', sep='\t')
        print(ccle_df.dtypes)
        #sns.swarmplot(data=ccle_df, x='tissue', y='outliers')
        sns.boxplot(data=ccle_df,x='tissue', y='outliers', whis=np.inf, color=[]'white', 'red'])
        #sns.violinplot(data=ccle_df, x='tissue', y='outliers', inner=None)
        plt.xticks(rotation='vertical')
        plt.title('Outlier counts for each tissue type with dPSI > 10')
        plt.show()


def main(myCommandLine=None):
    '''
    Creates an object of the ourputAnalyzer class and passes command line arguments to that classs

    '''
    myFileReader = outputAnalyzer()
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
