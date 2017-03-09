import os
import scipy
import math
from math import *
import numpy
import numpy as np
from numpy import NaN, Inf, arange, isscalar, asarray, array
from scipy.io.wavfile import write,read
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import pyplot



# def find_enpoints(y):
# 	y = np.asarray(y,'float')
# 	y = y/max(y)
# 	window_len = 3501
# 	t = np.linspace(0,len(y)/Fs,len(y)+1)
# 	energy = 0;
# 	w = np.hamming(window_len)
# 	energy = np.convolve(y**2,w**2,'same')
# 	energy = energy/max(energy)
# 	thresh = 0.008;
# 	energy_thresh = (energy > thresh).astype('float')
# 	points = np.nonzero(abs(energy_thresh[1:] - energy_thresh[0:-1]))[0]
# 	start_points = [points[2*i] for i in range(len(points)/2)]
# 	end_points = [	points[2*i+1] for i in range(len(points)/2)]
# 	f,axarr = plt.subplots(3)
# 	axarr[0].plot(y)
# 	axarr[0].set_title('Original Signal')
# 	axarr[1].plot(energy)
# 	axarr[1].set_title('Windowed Energy')
# 	te = (energy_thresh[1:] - energy_thresh[0:-1])
# 	ce = np.zeros(len(te))
# 	pe = np.zeros(len(te))
# 	for i in range(len(te)):
# 		if(te[i] == 1):
# 			ce[i] = 1
# 		if(te[i] == -1):
# 			pe[i] = 1

# 	axarr[2].plot(ce)
# 	axarr[2].hold
# 	axarr[2].plot(pe)
# 	axarr[2].set_title('Endpoints')
# 	plt.show()



# 	return start_points,end_points


# female_sounds=os.listdir("../data/Digits female 8Khz Updated")
# name = female_sounds[0]
# Fs,y = read('../data/Digits female 8Khz Updated/'+name)
# y = y[2000:]
# start_points,end_points = find_enpoints(y)



Fs= 8000.0
frame_size= 0.025
frame_hop=0.010
lowerf= 300.0
upperf= Fs/2
no_of_filters=26
n_fft = 1024


def melFilter(no_of_filters,lowerf,upperf,Fs,no_of_DFT_points):
# // The function returns the mel filterbank gain function for positive frequency. 
# //Inputs:
# //no_filters : no. of filters
# //lowerf : lower bound on frequency(Hz) 
# //upperf : upper bound on frequency(Hz)
# //Fs : Sampling frequency(samples/s)
# //no_of_DFT_points : no. of DFT points used in FFT calculation
# //Output:
# //filter_bank: contains the freuency response of filterbanks.
	no_of_DFT_points = n_fft
	mel_freq_lower = 2595*log10(1+(lowerf/700));
	mel_freq_upper = 2595*log10(1+(upperf/700));
	mel_spacing = (mel_freq_upper - mel_freq_lower)/(no_of_filters+1);
	mel_filter_cutoffs = [mel_freq_lower + i*mel_spacing for i in range(no_of_filters + 2)];
	linear_filter_cutoffs = [(10**(mel_filter_cutoffs[i]/2595) - 1)*700 for i in range(len(mel_filter_cutoffs))]
	linear_filter_cutoffs_discrete = [np.round(i*(no_of_DFT_points)/Fs) for i in linear_filter_cutoffs];

	filter_bank = np.zeros([no_of_filters,no_of_DFT_points/2+1]); 
	for i in range(1,len(linear_filter_cutoffs_discrete)-1):
	    for j in range(int(linear_filter_cutoffs_discrete[i-1]),int(linear_filter_cutoffs_discrete[i+1])+1):
	        if (j == linear_filter_cutoffs_discrete[i]):
	            filter_bank[i-1,j] = 1

	        elif(j < linear_filter_cutoffs_discrete[i]):
	            filter_bank[i-1,j] = (j-linear_filter_cutoffs_discrete[i-1])/(linear_filter_cutoffs_discrete[i] - linear_filter_cutoffs_discrete[i-1]);
	        else:
	            filter_bank[i-1,j] = (linear_filter_cutoffs_discrete[i+1]-j)/(linear_filter_cutoffs_discrete[i+1] - linear_filter_cutoffs_discrete[i]);
	  

	for k in range(no_of_filters):
		plt.plot(filter_bank[k,:])
		plt.hold
		plt.draw()
		

	plt.show()
 	return filter_bank   

filter_bank = melFilter(no_of_filters,lowerf,upperf,Fs,n_fft)
