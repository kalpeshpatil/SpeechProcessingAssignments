
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



def find_enpoints(y):
	y = np.asarray(y,'float')
	y = y/max(y)
	window_len = 3501
	t = np.linspace(0,len(y)/Fs,len(y)+1)
	energy = 0;
	w = np.hamming(window_len)
	energy = np.convolve(y**2,w**2,'same')
	energy = energy/max(energy)
	thresh = 0.008;
	energy_thresh = (energy > thresh).astype('float')
	points = np.nonzero(abs(energy_thresh[1:] - energy_thresh[0:-1]))[0]
	start_points = [points[2*i] for i in range(len(points)/2)]
	end_points = [	points[2*i+1] for i in range(len(points)/2)]
	return start_points,end_points

#pre_emphasis filter
def pre_emphasis_filter(x,alpha):
    y = np.zeros(len(x))
    y[0] = x[0]
    for i in range(1,len(x)):
        y[i]= x[i] - alpha*x[i-1]
    return y

male_sounds=os.listdir("../data/Digits male 8Khz Updated")
if (not os.path.isdir("Male_segmented")):
	os.mkdir('Male_segmented')

female_sounds=os.listdir("../data/Digits female 8Khz Updated")
if (not os.path.isdir("Female_segmented")):
	os.mkdir('Female_segmented')

def segment_digits():
	digit_list = ['zero','one','two','three','four','five','six','seven','eight','nine']
	for name in male_sounds:
		Fs,y = read('../data/Digits male 8Khz Updated/'+name)
		y = y[2000:]
		start_points,end_points = find_enpoints(y)

		for p in range(len(end_points)):
			digit =pre_emphasis_filter(y[start_points[p]:end_points[p]],0.95)
			digit = np.asarray(digit,'float')
			digit /= np.max(np.abs(digit))

			if (not os.path.isdir("Male_segmented/"+name[13:-5])):
					os.mkdir('Male_segmented/'+name[13:-5])
			write('Male_segmented/'+name[13:-5]+'/'+digit_list[int(floor(p/2))]+'_'+str(2*int(name[-5])-p%2)+'.wav',Fs,digit)



	for name in female_sounds:
		Fs,y = read('../data/Digits female 8Khz Updated/'+name)
		y = y[2000:]
		start_points,end_points = find_enpoints(y)

		for p in range(len(end_points)):
			digit =pre_emphasis_filter(y[start_points[p]:end_points[p]],0.95)
			digit = np.asarray(digit,'float')
			digit /= np.max(np.abs(digit))
			if (not os.path.isdir("Female_segmented/"+name[13:-5])):
					os.mkdir('Female_segmented/'+name[13:-5])
			write('Female_segmented/'+name[13:-5]+'/'+digit_list[int(floor(p/2))]+'_'+str(2*int(name[-5])-p%2)+'.wav',Fs,digit)

if __name__=='__main__':
    segment_digits()