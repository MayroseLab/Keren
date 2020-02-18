from scipy.stats import spearmanr, wilcoxon, ttest_rel
from sklearn.metrics import mean_squared_error
from Bio import AlignIO, SeqIO
from alterAlignments import get_format
import re
import numpy as np


#### measurement calculations #####

# auxiliary function to calculate the mean absolute value of error
def mrad(predictions, targets):
	if len(predictions) != len(targets):
		print("error from mrad calculation: tagets vector and prediction vector lengths differ")
		return "NA"
	mrad_val = 0
	pos_num = len(predictions)
	for i in range(pos_num):
		try:
			mrad_val += abs(((predictions[i] - targets[i])/targets[i]))
		except:
			continue
	mrad_val = mrad_val / pos_num
	return mrad_val

# auxiliary function to calculate the mean relative mse
def rmse(predictions, targets):
	if len(predictions) != len(targets):
		print("error from rmse calculation: tagets vector and prediction vector lengths differ")
		return "NA"
	rmse = 0
	pos_num = len(predictions)
	for i in range(pos_num):
		try:
			rmse += (((predictions[i] - targets[i])/targets[i])**2)
		except:
			continue
	rmse = (rmse / pos_num)**0.5
	return rmse

# auxiliary function to calculate the mean mse
def mse(predictions, targets):
	if len(predictions) != len(targets):
		print("error from mse calculation: tagets vector and prediction vector lengths differ")
		return "NA"
	mse = 0
	pos_num = len(predictions)
	for i in range(pos_num):
		mse += (((predictions[i] - targets[i]))**2)
	mse = (mse / pos_num)**0.5
	return mse

# auxiliary function to calculate the ratio between the spaces and characters in the alignment
def get_spaces_to_chars_ratio(msa_path, ref_seq_name):
	ref_seq_catcher = re.compile(">" + ref_seq_name + "s+([^\n]+)")
	with open(msa_path, "r") as msa_file:
		msa_content = msa_file.read()
		regex_cought = ref_seq_catcher.search(msa_content)
		if regex_cought != None:
			orig_seq = ref_seq_catcher.search(msa_content).group(1)
		else:
			return "NA"
		adapted_seq = orig_seq.replace("\n", "")
		alignment_length = len(adapted_seq)
		num_of_sequences = msa_content.count(">")
		spaces_num = msa_content.count("-")
		chars_num = int(alignment_length*num_of_sequences)-spaces_num
	return float(spaces_num/float(chars_num))

# auxiliary function to calculate the spearman correlation between two lists
def spearman(predictions, targets):
	pred_vector = np.array(predictions)
	tar_vector = np.array(targets)
	return spearmanr(pred_vector, tar_vector)

# auxiliary function to calculate the wilcoxon correlation coefficient between two vectors
def wilcoxon(vector_1, vector_2):
	T, p_val = wilcoxon(np.array(vector_1), np.array(vector_2))
	return p_val

# auxiliary function to calculate the non relative absolute error
def abs_dev(predictions, targets):
	if len(predictions) != len(targets):
		print("error from mse calculation: tagets vector and prediction vector lengths differ")
		return "NA"
	abs_err = 0
	pos_num = len(predictions)
	for i in range(pos_num):
		abs_err += (((predictions[i] - targets[i]))**2)
	abs_err = (abs_err / pos_num)**0.5
	return abs_err
