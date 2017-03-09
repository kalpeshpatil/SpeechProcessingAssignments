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
from mfcc_features import *
from scipy import spatial
import pickle

fs=8000
d=['zero','one','two','three','four','five','six','seven','eight','nine']
male_sounds=os.listdir("../data/Digits male 8Khz Updated")


male_speakers=os.listdir("Male_segmented")
female_speakers=os.listdir("Female_segmented")


all_speakers = male_speakers + female_speakers

CodeBook = {}
n_frame_dict = {}

for digit in d:
	CodeBook[digit] = {}
	n_frame_dict[digit] = {}

	for speaker in all_speakers:
		if(speaker in male_speakers):
			fp = "Male_segmented/"+speaker
		if(speaker in female_speakers):
			fp = "Female_segmented/"+speaker

		CodeBook[digit][speaker] = []
		n_frame_dict[digit][speaker] = []
		for tt in range(1,5):		
			Fs,signal = read(fp+'/'+digit+'_'+str(tt)+'.wav')
			feature_mat = extract_features(signal)
			n_frame_dict[digit][speaker].append(feature_mat.shape[1])
			for i in range(feature_mat.shape[1]):
				CodeBook[digit][speaker].append(feature_mat[:,i])

	print digit


# pickle.dump(CodeBook,open('CodeBook_v1','wb'))
# pickle.dump(n_frame_dict,open('n_frame_dict_v1','wb'))

CodeBook = pickle.load(open('CodeBook_v1','r'))
n_frame_dict = pickle.load(open('n_frame_dict_v1','r'))

confusion_matrix = np.zeros([10,10])

# test_speaker = male_speakers[0]

for test_speaker in all_speakers:
	train_speakers = male_speakers+female_speakers
	train_speakers.remove(test_speaker)

	for test_digit in d:
		for utterance in range(1,5):

			n_frames = n_frame_dict[test_digit][test_speaker]
			test_mat = CodeBook[test_digit][test_speaker][sum(n_frames[0:(utterance-1)]):sum(n_frames[0:utterance])]

			sum_dist = np.zeros(10)
			for digit in d:
				for l in range(len(test_mat)):
					test_vec = np.asarray(test_mat)[l,:]	
					min_dist = float('Inf')
					for speaker in train_speakers:
						temp_list = CodeBook[digit][speaker]
						curr_dist,index = spatial.KDTree(temp_list).query(test_vec)
						if(curr_dist < min_dist):
							min_dist = curr_dist

					sum_dist[d.index(digit)]+=min_dist

			pred_digit = np.argmin(sum_dist)
			print test_speaker,pred_digit,test_digit
			confusion_matrix[d.index(test_digit),pred_digit] += 1

	np.save('BOF_confusion_matrix',confusion_matrix)


wer = 1 - np.trace(confusion_matrix)/640
print wer