#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os
import time
import argparse
import pandas as pd
import json
import get_filename
import collections
from collections import defaultdict

def write_log(locus, log_text):
	logfile = get_filename.logfile(locus)
	with open(logfile, 'a') as f:
		f.write(log_text)

def options():
	parser = argparse.ArgumentParser(description='PhyloCode conversion')
	parser.add_argument(
		'-s', '--snp',
		type = str,
		dest = 'snp_file',
		required = True, 
		help = 'SNP distance file'
	) 
	parser.add_argument(
		'-c', '--circle',
		type = str,
		dest = 'circle_list',
		required = True, 
		help = 'circle file'
	) 
	parser.add_argument(
		'-o', '--out',
		type = str,
		dest = 'out_dir',
		default = '',
		help = 'Out put directory name'
	) 

	args = parser.parse_args()
	sdis_out = args.snp_file
	circle_list = args.circle_list
	out_dir = args. out_dir
	return(sdis_out, circle_list, out_dir)

class UnionFind:
	def __init__(self):
		self.parent = {}
		
	def find(self, x):
		if self.parent.setdefault(x, x) != x:
			self.parent[x] = self.find(self.parent[x])
		return self.parent[x]
	
	def union(self, x, y):
		rootX = self.find(x)
		rootY = self.find(y)
		if rootX != rootY:
			self.parent[rootY] = rootX

def group_and_merge_arrays(arrays_dic):
	uf = UnionFind()
	
	# Make a dictionary, its keys are the numbers of the list and the values are indices of the list.
	number_to_index = {}
	
	# Numbers in lists are joined into same group using UnionFind.
	for i, array in arrays_dic.items():
		for num in array:
			if num in number_to_index:
				uf.union(i, number_to_index[num])
			number_to_index[num] = i

	# Lists are joined based on groups each index is belonging to. 
	groups = defaultdict(list)
	for i, array in arrays_dic.items():
		root = uf.find(i)
		groups[root].extend(array)

	# Return results without duplications
	return [sorted(set(group)) for group in groups.values()]

def snp_grouping(locus):
	def parse_allele_number(allele_id):
		try:
			return int(str(allele_id).rsplit("_", 1)[1])
		except Exception:
			raise ValueError(
				f"Allele ID must end with '_<number>', but got: {allele_id}"
			)

	write_log(locus, f"snp_grouping\tBegin\t{str(time.ctime())}\n")
	sdis_out = get_filename.sdis_out(locus)
	# read TSV file as pandas dataframe
	df_Sdis = pd.read_csv(sdis_out, index_col=0, sep='\t')
	key_list = df_Sdis.index.to_list()  #key_list = allele name list
	
	# Alleles showing distance of 1 are listed for each allele
	SNPdis_group_dic = {}
	closest_list = []
	for index, item in df_Sdis.iterrows():
		for key in key_list:
			if item[key] <= 1:
				closest_list.append(parse_allele_number(key))
		SNPdis_group_dic[parse_allele_number(index)] = closest_list.copy()
		closest_list.clear()
	
	group_list = group_and_merge_arrays(SNPdis_group_dic)
	
	selection_dic = {} #{allele : [allele group]}
	
	# Select most frequently observed allele 2025/5/26
	allele_count = []
	for group in group_list:
		allele_count.clear()
		for allele in group:
			allele_count.extend(SNPdis_group_dic[allele])
		represent_allele = collections.Counter(allele_count).most_common()[0][0]
		selection_dic[f"{locus}_{represent_allele}"] = [f"{locus}_{name}" for name in group]
	

	dir_path = get_filename.dir_path(locus)
	with open(f"{dir_path}/{locus}_selection.json", 'w') as f:
		json.dump(selection_dic, f)
	write_log(locus, f"snp_grouping\tEnd\t{str(time.ctime())}\n")
	return(selection_dic, df_Sdis, key_list)

def snpgroup_list(locus, circle_list, out_dir, selection, df_Sdis, key_list):
	with open(circle_list, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]
	
	# Combine neighboring three groups
	trigroup = []
	trigroup_dict = {}
	trigroup_all = []
	for n in range(len(circle_list)):
		trigroup.append(circle_list[n-2])
		trigroup.append(circle_list[n])
		trigroup.extend(selection[circle_list[n-1]])
		trigroup_dict[circle_list[n-1]] = trigroup.copy()
		trigroup.clear()
		trigroup_all.append(circle_list[n-2])
		trigroup_all.append(circle_list[n])
		trigroup_all.extend(selection[circle_list[n-1]])
	
	# Check and append remaining alleles
	remain_list = []
	for allele in key_list:
		if allele in trigroup_all:
			pass
		else:
			remain_list.append(allele)
	if remain_list:
		df_circle = df_Sdis.loc[remain_list,circle_list]
		for allele in remain_list:
			trigroup_dict[df_circle.loc[allele].idxmin()].append(allele)
	
	#Remove group with 3 alleles
	keys_to_remove = []
	for ky, val in trigroup_dict.items():
		if len(val) <= 3:
			keys_to_remove.append(ky)

	for key in keys_to_remove:
		del trigroup_dict[key]

	#Save to list file
	for key, value in trigroup_dict.items():
		file_name = os.path.join(out_dir, key + ".list")
		with open(file_name, "w") as file:
			file.write('\n' .join(value))
