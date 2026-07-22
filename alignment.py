#!/usr/bin/python3
# -*- coding:utf-8 -*-

import time
import csv
from Bio import SeqIO
import pandas as pd
import copy
import os
import subprocess
import get_filename

def write_log(locus, log_text):
	logfile = get_filename.logfile(locus)
	with open(logfile, 'a') as f:
		f.write(log_text)

def fasta_select(in_fasta, out_fasta, seq_names):
	seq_list = []
	if seq_names:
		with open(seq_names, 'r') as sn:
			seqs = sn.read()
			for each in seqs.split('\n'):
				seq_list.append(each)

	with open(in_fasta, 'r') as ifas, open(out_fasta, 'w') as ofas:
		seq = SeqIO.parse(ifas, "fasta")
		for each in seq:
			if each.id in seq_list:
				SeqIO.write(each, ofas, 'fasta')

def blast_search_properly_seq(in_fasta, locus, blast_db, identity):
	# Open and read Fasta
	with open(in_fasta, "r") as fasta_file:
		records = list(SeqIO.parse(fasta_file, "fasta"))

		# Confirm the file is not empty
		if not records:
			print("Warning! Fasta file is empty")
			return

		# Get 1st sequence as representative
		first_sequence = records[0].seq
		sequence_length = len(first_sequence)
		
		#current_directory = os.getcwd()
		script_path = os.path.abspath(__file__)
		script_directory = os.path.dirname(script_path)
		database_path = blast_db
		bst_file = get_filename.bst_file(locus) #f"{locus}/{locus}_blast.tsv"
		output_file_path = bst_file

		# select word size based on the sequence length of 100 'modified from 57 to 100
		if sequence_length >= 100:
			blast_cmd = [
				"blastn",
				"-db", database_path,
				"-query", in_fasta,
				"-out", output_file_path,
				"-qcov_hsp_perc", "95",
				"-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
				"-num_alignments", "10"
			]
		else:
			blast_cmd = [
				"blastn",
				"-db", database_path,
				"-query", in_fasta,
				"-out", output_file_path,
				"-qcov_hsp_perc", "95",
				"-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
				"-num_alignments", "10",
				"-word_size", "11"
			]
		subprocess.run(blast_cmd, check=True, text=True)

def get_best_hits_by_query(blast_rows):
	best = {}
	for row in blast_rows:
		if len(row) < 11:
			continue

		query_id = row[0]
		bitscore = float(row[10])

		if query_id not in best or bitscore > float(best[query_id][10]):
			best[query_id] = row

	return best

def hit_select(locus, identity, in_fasta):
	directory_path = get_filename.dir_path(locus)
	bst_file = get_filename.bst_file(locus)
	hit_tsv = get_filename.hit_tsv(locus)
	hit_list_file = get_filename.hit_list(locus)

	nohit_list_file = os.path.join(directory_path, f"{locus}_nohit.list")
	excluded_tsv_file = os.path.join(directory_path, f"{locus}_excluded_blast.tsv")

	# All allele IDs in the original FASTA
	all_query_ids = []
	with open(in_fasta, "r") as fasta_handle:
		for record in SeqIO.parse(fasta_handle, "fasta"):
			all_query_ids.append(record.id)

	all_query_set = set(all_query_ids)

	# Read BLAST result TSV
	with open(bst_file) as fb:
		reader = csv.reader(fb, delimiter='\t')
		blast = [row for row in reader]

	hit_list = []
	Species_hitlist = []
	excluded_rows = []

	best_hits = get_best_hits_by_query(blast)
	seen_query = set(best_hits.keys())

	for query_id, row in best_hits.items():

		query_id = row[0]
		subject_id = row[1]
		pident = float(row[2])
		align_len = float(row[3])
		qlen = float(row[4])
		qcov = 100.0 * align_len / qlen if qlen else 0.0
		bitscore = float(row[10])

		if pident >= identity:
			hit_list.append([query_id, subject_id, pident])
			Species_hitlist.append(query_id)
		else:
			excluded_rows.append([
				query_id,
				subject_id,
				pident,
				qcov,
				bitscore,
				f"identity_below_{identity}"
			])

	# Queries with no BLAST output at all
	no_blast_hit = sorted(all_query_set - seen_query)

	for query_id in no_blast_hit:
		excluded_rows.append([
			query_id,
			"No_Hit",
			0,
			0,
			0,
			"no_blast_hit"
		])

	print(*Species_hitlist, sep='\n')

	# Accepted hit CSV
	col = ['', 'TOP_HIT', 'Sequence_Identity']
	df = pd.DataFrame(hit_list, columns=col)
	df = df.set_index('')
	df.to_csv(hit_tsv, sep='\t')

	# Accepted allele list
	with open(hit_list_file, 'w') as f:
		for hit_locus in Species_hitlist:
			f.write(f"{hit_locus}\n")

	# Detailed excluded table
	with open(excluded_tsv_file, 'w') as f:
		#f.write("allele\tTOP_HIT\tpident\treason\n")
		f.write("allele\tTOP_HIT\tpident\tqcov\tbitscore\treason\n")
		for row in excluded_rows:
			f.write("\t".join(map(str, row)) + "\n")

	# Simple no-hit/excluded allele list
	with open(nohit_list_file, 'w') as f:
		for row in excluded_rows:
			f.write(f"{row[0]}\n")

	write_log(
		locus,
		f"hit_select\taccepted={len(Species_hitlist)}\texcluded={len(excluded_rows)}\t"
		f"no_blast_hit={len(no_blast_hit)}\t{str(time.ctime())}\n"
		)

def validate_retained_alleles_after_blast(locus, min_alleles=3):
	hit_list_file = get_filename.hit_list(locus)

	with open(hit_list_file, "r") as f:
		retained = [line.strip() for line in f if line.strip()]

	if len(retained) < min_alleles:
		message = (
			f"Only {len(retained)} allele(s) were retained after BLAST filtering. "
			f"At least {min_alleles} retained alleles are required for PhyloCode construction."
		)
		write_log(locus, f"validate_retained_alleles_after_blast\tFailed\t{message}\t{str(time.ctime())}\n")
		raise ValueError(message)

	write_log(
		locus,
		f"validate_retained_alleles_after_blast\tEnd\t"
		f"retained={len(retained)}\t{str(time.ctime())}\n"
	)

def chk_loci(in_fasta, locus, blast_db, identity):
	write_log(locus, f"chk_loci\tBegin\t{str(time.ctime())}\n")
	print('Select loci belonging to target species.')
	blast_search_properly_seq(in_fasta, locus, blast_db, identity)
	hit_select(locus, identity, in_fasta)
	validate_retained_alleles_after_blast(locus, min_alleles=3)
	selected_fasta = get_filename.selected_fasta(locus)
	hit_list = get_filename.hit_list(locus)
	fasta_select(in_fasta, selected_fasta, hit_list)
	write_log(locus, f"chk_loci\tEnd\t{str(time.ctime())}\n")

def multiple_alignment(locus):
	write_log(locus, f"alignment\tBegin\t{str(time.ctime())}\n")
	selected_align = get_filename.selected_align(locus)
	selected_fasta = get_filename.selected_fasta(locus)
	hit_list = get_filename.hit_list(locus)
	with open(selected_align, "w") as out_handle:
		subprocess.run(
			["mafft", "--auto", selected_fasta],
			stdout=out_handle,
			check=True,
			text=True
		)
	write_log(locus, f"alignment\tEnd\t{str(time.ctime())}\n")

def count_lines(file_path):
	with open(file_path, 'r') as file:
		lines = file.readlines()
		return(len(lines))

def chk_snpfile(sdis_out, num_seq):
	# read CSV file as pandas dataframe
	df_Sdis = pd.read_csv(sdis_out, index_col=0, sep='\t')
	key_list = df_Sdis.index.to_list()  #key_list = allele name list
	if len(key_list) != num_seq:
		return(-1)
	else:
		return(0)

def snp_distance(locus):
	write_log(locus, f"snp_distance\tBegin\t{str(time.ctime())}\n")
	hit_list = get_filename.hit_list(locus)
	num_seq = count_lines(hit_list)
	selected_align = get_filename.selected_align(locus)
	sdis_out = get_filename.sdis_out(locus)
	ret_val = -1
	while ret_val == -1: # Avoid "snp-dists" ERROR by comparing the number of alleles in sdis_out and num_seq
		with open(sdis_out, "w") as out_handle:
			subprocess.run(
				["snp-dists", selected_align],
				stdout=out_handle,
				check=True,
				text=True
			)
		ret_val = chk_snpfile(sdis_out, num_seq)
	write_log(locus, f"snp_distance\tEnd\t{str(time.ctime())}\n")
