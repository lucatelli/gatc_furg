#!/usr/bin/env python
"""
Module with utility functions for the GATC @ FURG.

This module brings together a set of utility functions for astrophysics
that are developed as a subproduct of the research made by the Grupo de
Astrofísica Teórica e Computacional @ FURG. All functions presented here
are not directly related to each other and are developed by several authors
within the research group.

"""

import os
import numpy as np
import pyfits as pf

__author__ = "Fabricio Ferrari, Geferson Lucatelli and Leonardo Ferreira"
__credits__ = ["GATC Group"]
__version__ = "0.1"
__maintainer__ = "GATC"
__email__ = "leonardo.ferreira@furg.br"
__status__ = "Development"

def normalizePSFs(path):
	"""
		Normalize all PSFs in given path.

		If, for some reason, you got a folder with a bunch of unormalized
		PSFs, this function will normalize them. 
		
	"""
	psfsNames = sorted(os.listdir(path))

	for name in psfsNames:
		print name
		try:
			psf, hdr = pf.getdata(path + name, header=1)
			psf = psf / psf.sum()
			pf.writeto(path + name, psf, hdr, clobber=1)
		except:
			print "Unexpected error:", sys.exc_info()[0]																																																																																																																																																																																																																																																																																																																																																																																																																																																											
			continue
	
