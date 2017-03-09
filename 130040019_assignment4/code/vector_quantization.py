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

fs=8000

d=['zero','one','two','three','four','five','six','seven','eight','nine']
male_sounds=os.listdir("../data/Digits male 8Khz Updated")


male_speakers=os.listdir("Male_segmented")
female_speakers=os.listdir("Female_segmented")


all_speakers = male_speakers + female_speakers

CodeBook = {}
n_frame_dict = {}

CodeBook = pickle.load(open('CodeBook_v1','r'))
n_frame_dict = pickle.load(open('n_frame_dict_v1','r'))

n_clusters = 16


VQCodeBook = {}
centroids = {}





confusion_matrix = np.zeros([10,10])

# test_speaker = male_speakers[0]

for test_speaker in all_speakers:
	train_speakers = male_speakers+female_speakers
	train_speakers.remove(test_speaker)

	for digit in d:
		VQCodeBook[digit] = []
		centroids[digit] = []
		for speaker in train_speakers:
			mat = np.asarray(CodeBook[digit][speaker])
			if(VQCodeBook[digit] == []):
				VQCodeBook[digit] = mat
			else:
				VQCodeBook[digit] = np.concatenate([VQCodeBook[digit],mat],0)

		centroids[digit] = (kmeans(VQCodeBook[digit],n_clusters))[0]

	for test_digit in d:
		for utterance in range(1,5):

			n_frames = n_frame_dict[test_digit][test_speaker]
			test_mat = CodeBook[test_digit][test_speaker][sum(n_frames[0:(utterance-1)]):sum(n_frames[0:utterance])]

			sum_dist = np.zeros(10)
			for digit in d:
				for l in range(len(test_mat)):
					test_vec = np.asarray(test_mat)[l,:]	
					temp_list = np.asarray(centroids[digit])
					min_dist,index = spatial.KDTree(temp_list).query(test_vec)
					sum_dist[d.index(digit)]+=min_dist

			pred_digit = np.argmin(sum_dist)
			print test_speaker,pred_digit,test_digit
			confusion_matrix[d.index(test_digit),pred_digit] += 1.0


	np.save('VQ_confusion_matrix_n_cluster'+str(n_clusters),confusion_matrix)

# print confusion_matrix/64
wer = 1 - np.trace(confusion_matrix)/640

confusion_matrix = confusion_matrix/64
d = np.around(confusion_matrix,decimals = 3)
np.savetxt("../report/VQ_confusion_matrix"+str(n_clusters)+".csv", d, fmt = '%s')

print wer