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
from scipy.cluster.vq import vq,kmeans
import pickle

def find_dtw_distance(test_pattern,ref_pattern):
	n = test_pattern.shape[1];
	m = ref_pattern.shape[1];
	distMat = np.zeros([n, m]);

	for i in range(n):
		for j in range(m):

			distMat[i, j] = np.linalg.norm(np.subtract(test_pattern[:, i], ref_pattern[:, j]));

	DTW = np.zeros([n+1, m+1]);

	for i in range(1,n+1):
		DTW[i, 0] = float('Inf')

	for i in range(1,m+1):
		DTW[0, i] = float('Inf')

	DTW[0,0] = 0;

	for i in range(1,n+1):
		for j in range(1,m+1):
			cost = distMat[i-1, j-1];
			DTW[i, j] = cost + np.min([DTW[i-1, j], np.min([DTW[i-1, j-1], DTW[i, j-1]])]);

	return DTW[n, m];




fs=8000

d=['zero','one','two','three','four','five','six','seven','eight','nine']
male_sounds=os.listdir("../data/Digits male 8Khz Updated")


male_speakers=os.listdir("Male_segmented")
female_speakers=os.listdir("Female_segmented")


all_speakers = male_speakers + female_speakers
confusion_matrix = np.zeros([10,10])

CodeBook = {}
n_frame_dict = {}

CodeBook = pickle.load(open('CodeBook_v1','r'))
n_frame_dict = pickle.load(open('n_frame_dict_v1','r'))

for test_speaker in all_speakers:
	train_speakers = male_speakers+female_speakers
	train_speakers.remove(test_speaker)

	for test_digit in d:
		for utterance in range(1,5):

			n_frames = n_frame_dict[test_digit][test_speaker]
			test_mat = CodeBook[test_digit][test_speaker][sum(n_frames[0:(utterance-1)]):sum(n_frames[0:utterance])]

			sum_dist = np.zeros(10)
			min_dist = float('Inf')
			for digit in d:

				for speaker in train_speakers:
					for ref_utterance in range(1,5):
						n_ref_frames = n_frame_dict[digit][speaker]
						ref_mat = CodeBook[digit][speaker][sum(n_ref_frames[0:(ref_utterance-1)]):sum(n_ref_frames[0:ref_utterance])]

						test_pattern = np.transpose(np.asarray(test_mat))
						ref_pattern = np.transpose(np.asarray(ref_mat))
						curr_dist = find_dtw_distance(test_pattern,ref_pattern)
						if(curr_dist < min_dist):
							min_dist = curr_dist
							pred_digit = d.index(digit)

			print test_speaker,pred_digit,test_digit
			confusion_matrix[d.index(test_digit),pred_digit] += 1

	np.save('DTW_confusion_matrix',confusion_matrix)
# c = calcDTWDist(np.asarray(test_mat1),np.asarray(test_mat2))
wer = 1 - np.trace(confusion_matrix)/640

confusion_matrix = confusion_matrix/64
d = np.around(confusion_matrix,decimals = 3)
np.savetxt("../report/DTW_confusion_matrix.csv", d, fmt = '%s')

