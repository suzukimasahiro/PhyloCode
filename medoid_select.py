#!/usr/bin/python3
# -*- coding:utf-8 -*-

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def medoid_selection(dm_file, ls_file, max_taxa=80, linkage_method="average"):
	mt = np.loadtxt(dm_file, delimiter=",")

	with open(ls_file, "r") as f:
		label_list = [line.strip() for line in f if line.strip()]

	if mt.ndim != 2 or mt.shape[0] != mt.shape[1]:
		raise ValueError(f"Distance matrix must be square: {dm_file}")

	if mt.shape[0] != len(label_list):
		raise ValueError(
			f"Distance matrix size ({mt.shape[0]}) does not match "
			f"number of labels ({len(label_list)}) in {ls_file}"
		)

	if np.isnan(mt).any():
		raise ValueError(f"Distance matrix contains NaN: {dm_file}")

	if np.any(mt < 0):
		raise ValueError(f"Distance matrix contains negative values: {dm_file}")

	if not np.allclose(mt, mt.T, atol=1e-10):
		raise ValueError(f"Distance matrix is not symmetric: {dm_file}")

	np.fill_diagonal(mt, 0.0)

	n = len(label_list)

	if max_taxa < 1:
		raise ValueError(f"max_taxa must be >= 1, but got {max_taxa}")

	if max_taxa > 80:
		raise ValueError(f"max_taxa must be <= 80, but got {max_taxa}")

	if n == 0:
		raise ValueError(f"No labels found in {ls_file}")

	if n <= max_taxa:
		assignment = {label: [label] for label in label_list}
		return label_list, label_list, list(range(n)), assignment

	# If all distances are zero, hierarchical clustering is not meaningful.
	# Select the first max_taxa alleles deterministically.
	if np.allclose(mt, 0.0):
		medoid_num_list = list(range(max_taxa))
		medoid_list = [label_list[i] for i in medoid_num_list]
		assignment = {m: [m] for m in medoid_list}
		for i, label in enumerate(label_list[max_taxa:], start=max_taxa):
			assignment[medoid_list[i % max_taxa]].append(label)
		return label_list, medoid_list, medoid_num_list, assignment

	condensed = squareform(mt, checks=False)
	Z = linkage(condensed, method=linkage_method)
	cluster_labels = fcluster(Z, t=max_taxa, criterion="maxclust")

	medoid_indices = []
	assignment = {}

	for cluster_id in sorted(set(cluster_labels)):
		idx = np.where(cluster_labels == cluster_id)[0]
		subD = mt[np.ix_(idx, idx)]

		# medoid = allele with minimum total within-cluster distance
		medoid_local = int(np.argmin(subD.sum(axis=1)))
		medoid_index = int(idx[medoid_local])

		medoid_indices.append(medoid_index)
		medoid_name = label_list[medoid_index]
		assignment[medoid_name] = [label_list[i] for i in idx]

	medoid_indices = sorted(medoid_indices)
	medoid_list = [label_list[i] for i in medoid_indices]
	assignment = {m: assignment[m] for m in medoid_list}

	return label_list, medoid_list, medoid_indices, assignment
