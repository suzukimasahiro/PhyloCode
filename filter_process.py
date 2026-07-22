#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os

def medoid_group_processing(circle_list_file, out_dir, assignment, strict=False):
	os.makedirs(out_dir, exist_ok=True)

	with open(circle_list_file, "r") as f:
		circle_list = [line.strip() for line in f if line.strip()]

	if len(circle_list) == 0:
		raise ValueError(f"circle_list is empty: {circle_list_file}")

	missing = [m for m in circle_list if m not in assignment]
	if missing and strict:
		raise ValueError(
			f"{len(missing)} medoids in circle list are missing from assignment. "
			f"Examples: {missing[:5]}"
		)

	for i, medoid in enumerate(circle_list):
		prev_medoid = circle_list[i - 1]
		next_medoid = circle_list[(i + 1) % len(circle_list)]
		members = assignment.get(medoid, [medoid])

		group = [prev_medoid, next_medoid] + members

		unique_group = []
		for x in group:
			if x not in unique_group:
				unique_group.append(x)

		# Very small groups do not help local ordering, but keeping them is harmless.
		file_name = os.path.join(out_dir, medoid + ".list")
		with open(file_name, "w") as file:
			file.write("\n".join(unique_group))
			
