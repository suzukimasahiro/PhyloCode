#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os

def dir_path(locus):
	return(os.path.abspath(locus))

def result_path(locus):
	return(os.getcwd() + "/DATABASE")

def bst_file(locus):
	return(dir_path(locus) + f"/{locus}_blast.tsv")

def hit_tsv(locus):
	return(dir_path(locus) + f"/{locus}_hit.tsv")

def hit_list(locus):
	return(dir_path(locus) + f"/{locus}_hit.list")

def selected_align(locus):
	return(dir_path(locus) + f"/{locus}_selected.aln")

def selected_fasta(locus):
	return(dir_path(locus) + f"/{locus}_selected.fasta")

def sdis_out(locus):
	return(dir_path(locus) + f"/{locus}_Sdis.tsv")

def ss_list(locus): # representative alleles list file name
	return(dir_path(locus) + f"/{locus}_represent.list")

def tree_dir(locus):
	return(dir_path(locus) + "/Init_tree")

def logfile(locus):
	return(dir_path(locus) + f"/{locus}.log")
