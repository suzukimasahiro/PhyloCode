#!/usr/bin/python3
# -*- coding:utf-8 -*-

##### Dependencies #####
# ncbi blast
# FastTree or raxml-ng
# snp-dists
# mafft
# yukimayuli-gmz/lpnet
# python libraries
#   bio-python
#   pandas
#   numpy
#   matplotlib
#   scipy
#   rpy2
########################

import argparse
import os
import sys
import time
import gc
import shutil
from Bio import SeqIO
import get_filename
import alignment
import quantification


def write_log(locus, log_text):
	logfile = get_filename.logfile(locus)
	with open(logfile, 'a') as f:
		f.write(log_text)

def get_directory_and_file_names(relative_path):
	last_slash_index = relative_path.rfind('/')
	if last_slash_index == -1:
		directory_name = ''
		file_name = relative_path
	else:
		directory_name = relative_path[:last_slash_index]
		file_name = relative_path[last_slash_index + 1:]
	return directory_name, file_name

def validate_input_fasta(in_fasta, locus):
	"""
	Validate input allele FASTA before starting the PhyloCode workflow.

	Requirements:
	1. The FASTA file must contain at least 3 allele sequences.
	2. Each sequence ID must follow the format <locus>_<allele_number>.
	3. The prefix before the final underscore must match the locus name.
	4. The allele_number must be an integer.
	5. Sequence IDs must be unique.
	6. Sequences must not be empty.

	The validation result is written to the locus-specific log file.
	"""

	write_log(locus, f"validate_input_fasta\tBegin\t{str(time.ctime())}\n")

	records = list(SeqIO.parse(in_fasta, "fasta"))

	error_messages = []

	if len(records) < 3:
		error_messages.append(
			f"Input FASTA must contain at least 3 allele sequences, "
			f"but {len(records)} sequence(s) were found: {in_fasta}"
		)

	seen_ids = set()
	invalid_ids = []
	wrong_locus_ids = []
	duplicate_ids = []
	empty_sequence_ids = []

	for record in records:
		seq_id = record.id

		# Check duplicate IDs
		if seq_id in seen_ids:
			duplicate_ids.append(seq_id)
		else:
			seen_ids.add(seq_id)

		# Check empty sequence
		if len(record.seq) == 0:
			empty_sequence_ids.append(seq_id)

		# Check <locus>_<number> format
		if "_" not in seq_id:
			invalid_ids.append(seq_id)
			continue

		prefix, allele_num = seq_id.rsplit("_", 1)

		if prefix != locus:
			wrong_locus_ids.append(seq_id)

		try:
			int(allele_num)
		except ValueError:
			invalid_ids.append(seq_id)

	if duplicate_ids:
		error_messages.append(
			"Duplicate sequence IDs were found. "
			f"Examples: {duplicate_ids[:5]}"
		)

	if empty_sequence_ids:
		error_messages.append(
			"Empty sequences were found. "
			f"Examples: {empty_sequence_ids[:5]}"
		)

	if invalid_ids:
		error_messages.append(
			"Sequence IDs must follow the format <locus>_<allele_number>, "
			"where <allele_number> is an integer. "
			f"Examples of invalid IDs: {invalid_ids[:5]}"
		)

	if wrong_locus_ids:
		error_messages.append(
			f"Sequence ID prefixes must match the FASTA filename-derived locus name '{locus}'. "
			f"Examples of mismatched IDs: {wrong_locus_ids[:5]}"
		)

	if error_messages:
		log_text = (
			"validate_input_fasta\tFailed\t"
			f"{str(time.ctime())}\n"
			+ "\n".join(f"validate_input_fasta\tError\t{msg}" for msg in error_messages)
			+ "\n"
		)
		write_log(locus, log_text)

		raise ValueError(
			"Input FASTA validation failed:\n"
			+ "\n".join(f"- {msg}" for msg in error_messages)
		)

	success_message = (
		f"Input FASTA validation passed: "
		f"{len(records)} allele sequences were found for locus '{locus}'."
	)

	write_log(locus, f"validate_input_fasta\tEnd\t{success_message}\t{str(time.ctime())}\n")
	print(success_message)

def options():
	parser = argparse.ArgumentParser(description='PhyloCode conversion')
	parser.add_argument(
		'-i', '--in',
		type = str,
		dest = 'infile',
		required = True, 
		help = 'multi fasta'
	) 
	parser.add_argument(
		'-d', '--db',
		type = str,
		dest = 'ref_genome_db',
		required = True, 
		help = 'blast database including complete genome sequences to select loci belonging to the species sensu stricto'
	) 
	parser.add_argument(
		'-o', '--out',
		type = str,
		dest = 'out_dir',
		default = '',
		help = 'Output directory name'
	) 
	parser.add_argument(
		'--identity',
		type = int,
		dest = 'seq_identity',
		required = False,
		default = 95,
		help = 'threshold of sequence identity to select sequences belonging to target species'
	) 
	parser.add_argument(
		'--force',
		dest = 'force_write',
		action="store_true",
		default = False,
		help = 'Over write on existing files'
	) 
	parser.add_argument(
		'--discard_temp',
		dest = 'discard_temporary_files',
		action="store_true",
		default = False,
		help = 'discard temporary files'
	) 
	parser.add_argument(
		'--method',
		type = str,
		dest = 'method',
		choices = ['FastTree', 'RAxML'],
		default = 'FastTree',
		help = 'Phylogeny method. FastTree or RAxML. Default: FastTree'
	) 
	parser.add_argument(
		'--threads',
		type = int,
		dest = 'num_threads',
		required = False,
		default = 4,
		help = 'The number of threads for RAxML-NG'
	) 
	parser.add_argument(
		"--max-lpnet-taxa",
		type=int,
		dest="max_lpnet_taxa",
		default=80,
		help="Maximum number of medoid alleles used as lpnet input. Default: 80."
	)
	args = parser.parse_args()
	in_fasta = args.infile 
	blast_db = args.ref_genome_db
	out_dir = args.out_dir
	identity = args.seq_identity
	force_write = args.force_write
	rm_temp = args.discard_temporary_files
	method = args.method
	threads = args.num_threads
	max_lpnet_taxa = args.max_lpnet_taxa
	
	in_fasta_abs = os.path.abspath(in_fasta)
	blast_db_dir, blast_db_name = get_directory_and_file_names(blast_db)
	blast_db_dir_abs = os.path.abspath(blast_db_dir)
	blast_db_abs = f"{blast_db_dir_abs}/{blast_db_name}"

	return(in_fasta_abs, blast_db_abs, identity, out_dir, force_write, rm_temp, method, threads, max_lpnet_taxa)

# Initialize Dir
def init_dir(locus, force_write):
	print('Initializing directories')
	if force_write and os.path.isdir(locus):
		print("Remove previous output directory.")
		shutil.rmtree(locus)

	try:
		os.makedirs(locus)
	except FileExistsError:
		pass
	try:
		os.makedirs("DATABASE")
	except FileExistsError:
		pass
	try:
		os.makedirs(get_filename.tree_dir(locus))
	except FileExistsError:
		pass

def phylocode(in_fasta, blast_db, identity, locus, rm_temp, method, threads, max_lpnet_taxa):
	logfile = get_filename.logfile(locus)
	log_dic = {'chk_loci':'NotYet', 'alignment':'NotYet','snp_distance':'NotYet','phylocode':'NotYet'}

	#init_dir(locus, force_write)

	if os.path.isfile(logfile):
		with open(logfile, 'r') as file:
			for line in file:
				parts = line.rstrip("\n").split("\t")
				if len(parts) < 2:
					continue
				key, value = parts[0], parts[1]
				if key in log_dic:
					log_dic[key] = value

	if log_dic['chk_loci'] == 'End':# and force_write == False:
		write_log(locus, f"chk_loci\tEnd\t{str(time.ctime())}\n")
	else:
		#write_log(locus, f"Target locus\t{locus}\n")
		alignment.chk_loci(in_fasta, locus, blast_db, identity)
			# Check species of each allele to remove sequences not belonging to target species

	if log_dic['alignment'] == 'End':# and force_write == False:
		write_log(locus, f"alignment\tEnd\t{str(time.ctime())}\n")
	else:
		alignment.multiple_alignment(locus)
		# multiple alignment

	if log_dic['snp_distance'] == 'End':# and force_write == False:
		write_log(locus, f"snp_distance\tEnd\t{str(time.ctime())}\n")
	else:
		alignment.snp_distance(locus)

	if log_dic['phylocode'] == 'End':# and force_write == False:
		write_log(locus, f"phylocode\tEnd\t{str(time.ctime())}\n")
	else:
		quantification.calc_phylocode(locus, rm_temp, method, threads, in_fasta, max_lpnet_taxa)
	gc.collect()

if __name__ == "__main__":
	in_fasta, blast_db, identity, out_dir, force_write, rm_temp, method, threads, max_lpnet_taxa = options()

	if not 1 <= max_lpnet_taxa <= 80:
		print("--max-lpnet-taxa must be between 1 and 80.", file=sys.stderr)
		sys.exit(1)
		
	file_name = os.path.basename(in_fasta)
	locus, _ = os.path.splitext(file_name)

	if out_dir != '':
		try:
			os.makedirs(out_dir)
		except FileExistsError:
			pass
		os.chdir(out_dir)

	init_dir(locus, force_write)

	write_log(locus, f"Target locus\t{locus}\n")

	try:
		validate_input_fasta(in_fasta, locus)
	except ValueError as e:
		print(str(e), file=sys.stderr)
		sys.exit(1)

	phylocode(
		in_fasta,
		blast_db,
		identity,
		locus,
		rm_temp,
		method,
		threads,
		max_lpnet_taxa
	)
