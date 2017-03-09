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


Fs= 8000.0
frame_size= 0.025
frame_hop=0.010
lowerf= 300.0
upperf= Fs/2
no_of_filters=26
n_fft = 1024



# mel_filter_bank = np.zeros(no_of_filters,n_fft/2)
# mel_points = zeros(no_of_filters)
# temp = (mel(upperf) - mel(lowerf))/(no_of_filters + 1)
# mel_points = [lowerf + i*mel_spacing for i in range(no_of_filters+2)]

def melFilter(no_of_filters,lowerf,upperf,Fs,no_of_DFT_points):

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
	  

	# for k in range(no_of_filters):
	# 	plt.plot(filter_bank[k,:])
	# 	plt.hold
	# 	plt.draw()
		
	# plt.show()
 	return filter_bank   
    
filter_bank = melFilter(no_of_filters,lowerf,upperf,Fs,n_fft)


def extract_features(signal):
	n_frames= 1+ math.floor((len(signal)-(frame_size*Fs))/(frame_hop*Fs))
	MFCC = []
	for i in range(int(n_frames)):

		frame= signal[i*(frame_hop*Fs):(frame_size*Fs)+i*(frame_hop*Fs)]

		w = np.hamming(len(frame))
		w_frame = [w[i]*frame[i] for i in range(len(frame))]
		temp = np.fft.fft(w_frame,n_fft)
		dft = temp[0:(n_fft)/2+1	]
		dft_mag = [(abs(t))**2 for t in dft]
		mel_coefs = np.zeros(no_of_filters)
		mel_coefs = [np.sum([dft_mag[j]*filter_bank[t,j] for j in range(len(dft_mag))]) for t in range(no_of_filters)]
		# for tt in range(no_of_filters):
		# 	for jj in range(len(dft_mag)):
		# 		mel_coefs[tt] += dft_mag[jj]*filter_bank[tt,jj]

		log_mel_coefs = [20*log10(abs(coef)) for coef in mel_coefs]
		mfcc_temp = scipy.fftpack.dct(log_mel_coefs)
		mfcc = mfcc_temp[1:14]
		#mfcc=(abs(numpy.fft.ifft(mel_coefs)))[:13]  
		MFCC.append(mfcc)

	MFCC = np.asarray(MFCC)
	delta_vecs = np.zeros(MFCC.shape)
	for i in range(1,MFCC.shape[1]-1):
	    delta_vecs[:,i] = np.subtract(MFCC[:,i+1] ,MFCC[:,i-1])/2
	           
	delta_delta_vecs = np.zeros(MFCC.shape)
	for i in range(1,MFCC.shape[1]-1):
	    delta_delta_vecs[:,i] = np.subtract(delta_vecs[:,i+1], delta_vecs[:,i-1])/2

	final_feature_vecs = np.transpose(np.concatenate([MFCC,delta_vecs,delta_delta_vecs],1))

	return final_feature_vecs

# plt.plot(np.transpose(final_feature_vecs))
# plt.show()




# if __name__=='__main__': 
    
#     lol=0    
