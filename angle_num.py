#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import copy
import math
import json
import get_filename

def make_order_list(circleList_file, dm_dir, dmlist_dir, mode):
	with open(circleList_file, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]

	# Get all list files in the directory
	file_names = [f for f in os.listdir(dmlist_dir) if f.endswith('.list')]
	file_contents_dict = {}
	suffix_to_remove = "_dm"

	disbw_dict = {}
	order_list = []
	order_dict = {}


	# Grouping with each file
	for file_name in file_names:
		with open(os.path.join(dmlist_dir, file_name), 'r') as file:
			file_contents = file.read().splitlines()
			base_name = os.path.splitext(f'{file_name}')[0]
			if base_name.endswith(suffix_to_remove):
				base_name = base_name[:-len(suffix_to_remove)]
			file_contents_dict[base_name] = file_contents
			csv_file_path = dm_dir + '/' + f'{base_name}_dm.csv'
			mt = np.loadtxt(csv_file_path, delimiter=',')
	
			s1 = circle_list.index(f'{base_name}')-1
			s2 = circle_list.index(f'{base_name}')+1
			if s2==len(circle_list):
				s2 = 0
			n1 = file_contents.index(circle_list[s1])
			n2 = file_contents.index(circle_list[s2])
			for seq in file_contents:
				num = file_contents.index(seq)
				dism = mt[num, n1]
				disp = mt[num, n2]
				disbw = dism - disp
				disbw_dict[seq] = disbw
			sorted_dict = {k: v for k, v in sorted(disbw_dict.items(), key=lambda item: item[1])}
			order_list = list(sorted_dict.keys())
			if len(order_list) >= 2:
				order_list.pop(0)
				order_list.pop(-1)
			else:
				order_list = []
			order_dict[base_name] = order_list.copy()
			disbw_dict.clear()
			order_list.clear()

	if mode == 'G2':
		order_list = []
		for allele in circle_list:
			#order_list.extend(order_dict[allele])
			if allele in order_dict:
				order_list.extend(order_dict[allele])
			else:
				order_list.append(allele)
		return(order_list)
	elif mode == 'G3G4':
		return(order_dict)
	else:
		order_list = []
		for allele in circle_list:
			if allele in order_dict.keys():
				order_list.extend(order_dict[allele])
			else:
				order_list.append(allele)
		return(order_list)

def allorder_listmk(circleList_file, dm_dir, dmlist_dir, orderList_file):
	allorder_list = make_order_list(circleList_file, dm_dir, dmlist_dir, 'allorder_listmk')
	with open(orderList_file, 'w') as file:
		for item in allorder_list:
			file.write(str(item) + '\n')

##### Update to fix order direction #####
def allele_number(allele_id):
	try:
		return int(str(allele_id).rsplit("_", 1)[1])
	except Exception:
		return str(allele_id)


def fix_orientation_by_allele_number(order):
	"""
	Fix orientation after anchoring the circular order.
	The first allele is kept fixed. The direction whose second allele has
	the smaller allele number is selected.
	"""
	if len(order) <= 2:
		return order

	forward = order
	reverse = [order[0]] + list(reversed(order[1:]))

	if allele_number(reverse[1]) < allele_number(forward[1]):
		return reverse
	return forward


def replace_with_cumulative_sum(dmcsv, label_list, order_list):
	dmt = np.loadtxt(dmcsv, delimiter=",")

	if len(order_list) == 0:
		raise ValueError("order_list is empty.")

	label_to_index = {label: i for i, label in enumerate(label_list)}

	missing = [x for x in order_list if x not in label_to_index]
	if missing:
		raise ValueError(
			f"{len(missing)} alleles in order_list are missing from label_list. "
			f"Examples: {missing[:5]}"
		)

	if len(order_list) == 1:
		return [0.0], order_list, [0.0], order_list[0]

	# Determine the largest adjacent gap in the circular order.
	bw_dis = []
	for v in range(len(order_list)):
		n1 = label_to_index[order_list[v]]
		n2 = label_to_index[order_list[v - 1]]
		bw_dis.append(float(dmt[n1, n2]))

	max_index = int(np.argmax(bw_dis))
	correct_order = order_list[max_index:] + order_list[:max_index]

	# Optional but recommended for reproducibility
	correct_order = fix_orientation_by_allele_number(correct_order)

	betw_dis = [0.0]
	for v in range(len(correct_order) - 1):
		m1 = label_to_index[correct_order[v]]
		m2 = label_to_index[correct_order[v + 1]]
		betw_dis.append(float(dmt[m1, m2]))

	result_list = []
	cumulative_sum = 0.0
	for dist in betw_dis:
		cumulative_sum += dist
		result_list.append(cumulative_sum)

	standard_allel = correct_order[0]
	return result_list, correct_order, betw_dis, standard_allel

def save_radian_dic(dis_dict, jsonfile):
	values = list(dis_dict.values())

	if len(values) == 0:
		raise ValueError("Distance dictionary is empty.")

	total = float(values[-1])

	if total == 0:
		radian_dict = {label: 0.0 for label in dis_dict.keys()}
	else:
		radian_dict = {
			label: float(value) * math.pi / total
			for label, value in dis_dict.items()
		}

	with open(jsonfile, "w") as file:
		json.dump(radian_dict, file, indent=2)

	return radian_dict

def draw_distr_map(radian_dict, locus, circleList_file):
	with open(circleList_file, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]

	# Draw a distribution map
	fig, ax = plt.subplots(figsize=(10, 6)) #(8, 4)
	ax.set_aspect('equal')

	# Plot semi-circle
	circle_x = np.cos(np.linspace(0, np.pi, 360))
	circle_y = np.sin(np.linspace(0, np.pi, 360))
	ax.plot(circle_x, circle_y, color='gray', linestyle='--')

	# Plot data on the semi-circle
	len_radian_dict = len (radian_dict)
	for label, angle in radian_dict.items():
		x = np.cos(angle)
		y = np.sin(angle)
		i = len_radian_dict % 6
		color_list = ['#6699ff', '#ffcc00', '#0099cc', '#ffcc99', '#99ccff', '#cc9933']
		ax.scatter(x, y, c=color_list[i])
		len_radian_dict += 1


	# Add data label
	for label in circle_list:
		angle = radian_dict[label]
		x = np.cos(angle)
		y = np.sin(angle)

		ax.scatter(x, y, c='red')

		offset_x = x / 6 
		offset_y = y / 6
		degree = math.degrees(angle)
		if degree >= 90:
			degree = degree + 180
		plt.text(x + offset_x, y + offset_y, f'{label} / {angle:.2f}', fontsize=8, ha='center', va='center', rotation=degree)

	# Define axis
	ax.set_xlim(-1.1, 1.1)
	ax.set_ylim(0, 1.1)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.legend()
	
	directory_path = get_filename.dir_path(locus)
	plt.savefig(f"{directory_path}/{locus}.eps", format='eps')

def G1_method(dmcsv, dmlist, circleList_file, jsonfile, locus):
	mt = np.loadtxt(dmcsv, delimiter=',')

	with open(dmlist, 'r') as f:
		label_list = [line.strip() for line in f if line.strip()]

	with open(circleList_file, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]
	
	standard_dis, correct_order, betw_dis, standard_allel = replace_with_cumulative_sum(dmcsv, label_list, circle_list)

	standard_dis_dict = dict(zip(correct_order, standard_dis))
	radian_dict = save_radian_dic(standard_dis_dict, jsonfile)
	draw_distr_map(radian_dict, locus, circleList_file)
	return(standard_allel)

def G2_method(dmcsv, dmlist, circleList_file, dm_dir, dmlist_dir, jsonfile, locus):
	dmt = np.loadtxt(dmcsv, delimiter=',')

	with open(dmlist, 'r') as f:
		label_list = [line.strip() for line in f if line.strip()]

	with open(circleList_file, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]

	order_list = make_order_list(circleList_file, dm_dir, dmlist_dir, 'G2')

	standard_dis, correct_order, betw_dis, standard_allel = replace_with_cumulative_sum(dmcsv, label_list, order_list)

	standard_dis_dict = dict(zip(correct_order, standard_dis))
	radian_dict = save_radian_dic(standard_dis_dict, jsonfile)
	draw_distr_map(radian_dict, locus, circleList_file)
	return(standard_allel)

def G3andG4_method(dmcsv, dmlist, circleList_file, dm_dir, dmlist_dir, jsonfile, locus, circleList):
	dmt = np.loadtxt(dmcsv, delimiter=',')

	with open(dmlist, 'r') as f:
		label_list = [line.strip() for line in f if line.strip()]

	with open(circleList_file, 'r') as f:
		circle_list = [line.strip() for line in f if line.strip()]

	order_dict = make_order_list(circleList_file, dm_dir, dmlist_dir, 'G3G4')

	
	sum_dis, correct_order, tmp_dis, tmp_standard = replace_with_cumulative_sum(dmcsv, label_list, circle_list)

	### Numbering to SNPgroup
	b_dis_dict = dict(zip(correct_order, sum_dis))

	## Split the list
	def split_list_around_element(input_list, delimiter):
		index = input_list.index(delimiter)

		# Make lists by splitting at the border value
		before_delimiter = input_list[:index]
		after_delimiter = input_list[index+1:]

		return before_delimiter, after_delimiter

	subb_dis_dict = {}

	for k, v in order_dict.items():
		dml_file_path = dmlist_dir + '/' + k + '_dm.list'
		with open(dml_file_path, 'r') as dml:
			sublabel_list = [line.strip() for line in dml if line.strip()]

		refnum = sublabel_list.index(k)

		csv_file_path = dm_dir + '/' + k + '_dm.csv'
		mt = np.loadtxt(csv_file_path, delimiter=',')
		m1 = circle_list.index(k)-1
		m2 = circle_list.index(k)+1
		if m2==len(circle_list):
			m2 = 0
		mnum1 = sublabel_list.index(circle_list[m1])
		mnum2 = sublabel_list.index(circle_list[m2])
		mmax = mt[mnum1, mnum2]
		Mnum1 = label_list.index(circle_list[m1])
		Mnum2 = label_list.index(circle_list[m2])
		Mmax = dmt[Mnum1, Mnum2]
		if mmax == 0:
			cor = 1.0
		else:
			cor = Mmax/mmax

		before, after = split_list_around_element(v, k)

		for mn in before:
			mndis = b_dis_dict[k] - mt[refnum, sublabel_list.index(mn)]*cor
			subb_dis_dict[mn] = mndis

		for pl in after:
			pldis = b_dis_dict[k] + mt[refnum, sublabel_list.index(pl)]*cor
			subb_dis_dict[pl] = pldis

	# Combine dictionary
	all_dis_dict = {**subb_dis_dict, **b_dis_dict}
	sorted_ad_dict = {k: v for k, v in sorted(all_dis_dict.items(), key=lambda item: item[1])}

	first_value = list(sorted_ad_dict.values())[0]
	sorted_ad_dict = {key: value - first_value for key, value in sorted_ad_dict.items()}

	standard_allel = list(sorted_ad_dict.keys())[0]

	radian_dict = save_radian_dic(sorted_ad_dict, jsonfile)
	draw_distr_map(radian_dict, locus, circleList)
	return(standard_allel)


