#!/usr/bin/python3
# -*- coding:utf-8 -*-

###yukimayuli-gmz/lpnet  required packages (Download the following packages to install lpnet.)

#CentOS
#harfbuzz-devel fribidi-devel libcurl-devel libxml2-devel freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel glpk-doc glpk-utils libglpk-dev

#Ubuntu
#libfontconfig1-dev libxml2-dev libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcurl4-openssl-dev glpk-doc glpk-utils libglpk-dev

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="rpy2.rinterface")

import numpy as np
import rpy2.rinterface_lib.callbacks
rpy2.rinterface_lib.callbacks.logger.setLevel("ERROR")
import rpy2.robjects as ro ### Requires Jinja2 3.x.x
import rpy2.interactive as r
import rpy2.interactive.packages
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
import medoid_select

rlib = importr("lpnet")

def lpnet(dm_file, ls_file, out_file):
	ro.r.assign("dm_file", dm_file)
	ro.r.assign("ls_file", ls_file)
	ro.r.assign("out_file", out_file)

	ro.r('M <- read.csv(dm_file, header =F)')
	ro.r('hoge <- matrix(as.matrix(M), nrow(M), ncol(M))')
	ro.r('taxaname <- c(scan(ls_file, what = character(), sep = "\n", blank.lines.skip = F))')
	ro.r('lpnet(hoge, tree.method = "nj", lp.package = "Rglpk", filename = out_file, taxaname = taxaname)')


def medoid_lpnet(dm_file, ls_file, out_file, max_taxa):
	label_list, medoid_list, medoid_num_list, assignment = medoid_select.medoid_selection(
		dm_file,
		ls_file,
		max_taxa=max_taxa
	)

	mt = np.loadtxt(dm_file, delimiter=",")
	np_array = mt[np.ix_(medoid_num_list, medoid_num_list)]

	with localconverter(default_converter + numpy2ri.converter):
		array = ro.conversion.py2rpy(np_array)

		ro.r.assign("sdismat", array)
		ro.r.assign("Rselection_list", medoid_list)
		ro.r.assign("out_file", out_file)

		ro.r('M <- as.matrix(sdismat)')
		ro.r('hoge <- matrix(as.matrix(M), nrow(M), ncol(M))')
		ro.r('taxaname <- c(Rselection_list)')
		ro.r('lpnet(hoge, tree.method = "nj", lp.package = "Rglpk", filename = out_file, taxaname = taxaname)')

	return label_list, medoid_list, assignment
