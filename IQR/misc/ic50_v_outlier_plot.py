
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.font_manager import FontProperties

def plot_correlation():
    colon = {'LoVo': [0.51, 9439], 'LS123': [0.42, 8752], 'SW837': [0.74, 12605]}
    blood = {'JeKo_1': [0.016, 9519], 'K_562': [0.39, 10227], 'RPMI_8226': [0.94, 8860], 'SR_786': [0.21, 13471]}

    plot(colon)
    plot(blood)

def plot(dict):
    a = []
    b = []
    for item in dict:
        a.append(dict[item][0])
        b.append(dict[item][1])

    print(a)
    print(b)
    plt.plot(a, b, c='blue', label=item)

    plt.autoscale(enable=True, axis='both', tight=None)
    plt.show()

def main():
    plot_correlation()

main()
