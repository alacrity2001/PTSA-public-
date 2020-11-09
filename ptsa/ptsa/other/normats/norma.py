# Authors: Paul Boniol, Michele Linardi, Federico Roncallo, Themis Palpanas
# Date: 08/07/2020
# copyright retained by the authors
# algorithms protected by patent application FR2003946
# code provided as is, and can be used only for research purposes




from .nA_normalmodel import *
from .nA_recurrent_sequences import *
from .tools import *
import matrix_profile as mp
import os
from ..matrix_profile.mp_wrapper import Join, selfJoin


class NormA():

	def __init__(self,pattern_length,nm_size):
		self.pattern_length = pattern_length
		self.nm_size = nm_size

	def run(self,ts_name,tot_length,
		percentage_sel=0.4,
		nm_name='tmp', clean = False):
		ts = get_bounded_ts(ts_name)[0:tot_length]

		if clean is True:
			ts = self.convert_zero(ts)
		else:
			pass
		recurrent_sequence,sequence_rec = extract_recurrent_sequences_random(ts, self.nm_size,percentage_sel=percentage_sel)
		listcluster,dendogram = clustering_method(recurrent_sequence)
		nms,scores_nms= choose_normalmodel(listcluster,recurrent_sequence, sequence_rec)
		
		self.normalmodel = [nms,scores_nms]

		if not os.path.exists(nm_name):
			os.makedirs(nm_name)

		for index_name,nm in enumerate(nms):
			with open(nm_name + '/' + str(index_name) ,"w") as f:
				for p in nm:
					f.write(str(p) + "\n")

		all_join = []
		for index_name in range(len(nms)):
			join,_ = Join(nm_name + '/' + str(index_name),ts_name,len(nms[index_name]),len(ts), self.pattern_length)
			join = np.array(join)
			if max(join) - min(join) == 0:
				pass
			else:
				join = (join - min(join))/(max(join) - min(join))
			all_join.append(join)
		
		join = [0]*len(all_join[0])
		for sub_join,scores_sub_join in zip(all_join,scores_nms):
			join = [float(j) + float(sub_j)*float(scores_sub_join) for j,sub_j in zip(list(join),list(sub_join))]

		join = np.array(join)
		join = running_mean(join,self.pattern_length)

		return join

	def run_motif(self,ts_name,tot_length,
		percentage_sel=0.4,
		nm_name='tmp'):
		ts = get_bounded_ts(ts_name)[0:tot_length]
		self_join,_ = selfJoin(ts_name, tot_length,self.pattern_length)
		recurrent_sequence,sequence_rec = extract_recurrent_sequences_motif(ts,self_join,self.nm_size,self.pattern_length)
		listcluster,dendogram = clustering_method(recurrent_sequence)
		nms,scores_nms = choose_normalmodel(listcluster,recurrent_sequence, sequence_rec)
		
		self.normalmodel = [nms,scores_nms]

		if not os.path.exists(nm_name):
			os.makedirs(nm_name)

		for index_name,nm in enumerate(nms):
			with open(nm_name + '/' + str(index_name) ,"w") as f:
				for p in nm:
					f.write(str(p) + "\n")

		all_join = []
		for index_name in range(len(nms)):
			join,_ = Join(nm_name + '/' + str(index_name),ts_name,len(nms[index_name]),len(ts), self.pattern_length)
			join = np.array(join)
			if max(join) - min(join) == 0:
				pass
			else:
				join = (join - min(join))/(max(join) - min(join))
			all_join.append(join)

		join = [0]*len(all_join[0])
		for sub_join,scores_sub_join in zip(all_join,scores_nms):
			join = [float(j) + float(sub_j)*float(scores_sub_join) for j,sub_j in zip(list(join),list(sub_join))]

		join = np.array(join)
		join = running_mean(join,self.pattern_length)

		return join
	
	def convert_zero(self,ts):
		data = ts.copy()
		for i in range(len(ts)):
			j = 0
			while data[i] == 0:
				if i+j < len(ts):
					data[i] = ts[i+j]
				else:
					data[i] = ts[i-1]
				j+= 1

		return data

