# !/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
## This is a main script that runs routines to create synthetic galaxies	##
## 						Created by Geferson Lucatelli 						##
##						lucatelli_asph_r2@virgophysics						##
##																			##
##		(original version 03.16.2016) (lastest version 04.06.2016)			##
##############################################################################

##############################################################################
##IMPORT SCRIPTS
##############################################################################
import galaxymakerlib
from galaxymakerlib import*
reload (galaxymakerlib)
import os

'''
Create synthetic galaxies, including in ellipticals and spirals components like bulges, disks and bars.
'''
'''
		The function galaxie() calls a specified profile, to create a synthetic galaxie. 
		returns an image array

		Some components that can be used to control the intensity of the profile in the galaxy are:
		ec		exponential contribution to the final galaxy intensity output
		buc		bulge contribution to the final galaxy intensity output
		dc 		disk contribution to the final galaxy intensity output
		barc	bar contribution to the final galaxy intensity output
		nsc		general sersic profile contribution to the final galaxy intensity output
		sct		tangent spiral contribution to the final galaxy intensity output
		scth	hyperbolic tangent spiral contribution to the final galaxy intensity output

		It can be created a simple synthetic galaxie in a more easy way simply calling a function.
'''
def create_galaxies_tan():
	'''
		This function join some functions of galaxymakerlib to create a more complex sythetic galaxie, 
		including spiral arms with the use of normal tangent function. 
		See paper http://arxiv.org/abs/0908.0892
		Returns an array image.
	'''
	nro__galaxias=1
	nsc=1.0
	In=30.0
	rn=40.0
	n=4.5
	ec=0.0
	Ie=35.0
	# re=np.linspace(35,70.0,nro__galaxias)
	re=30.0
	# sct=np.linspace(70,300.0,nro__galaxias)
	sct=270.0
	k=2.0
	p=1.0
	AA=1.5
	NN=3.3
	# NN=np.linspace(0.7,3.0,nro__galaxias)
	# phi0=np.linspace(0.1,3.0,nro__galaxias)
	phi0=2.9
	barc=0.0
	Ib=15.0
	rb=15.0
	nb=0.5
	qbar=0.2
	cbar=1.6
	q=0.95
	c=0.0
	i=0
	gal=spiral_galaxie(nsc,In,rn,n,ec,Ie,re,sct,k,p,AA,NN,phi0,barc,Ib,rb,nb,qbar,cbar,q,c)
	plot_save_image(gal,i,nsc,In,rn,n,ec,Ie,re,sct,k,p,AA,NN,phi0,barc,Ib,rb,nb,q,c)
	save_fits(gal,i)
	# for i in range(nro__galaxias):
	# 	gal=spiral_galaxie(nsc,In,rn,n,ec,Ie,re[i],sct[i],k,p,AA,NN[i],phi0[i],barc,Ib,rb,nb,q,c)
	# 	plot_save_image(gal,i,nsc,In,rn,n,ec,Ie,re[i],sct[i],k,p,AA,NN[i],phi0[i],barc,Ib,rb,nb,q,c)
	# 	save_fits(gal,i)
	# 	plt.clf()
create_galaxies_tan()
def call_mfmtk():
	data_fits=['gal_1.fits']
	for i in data_fits:
		os.system("python path "+str(i)+" ppsf.fits")
		print 'Gal ', str(i),' done.'
# call_mftk()
# def create_galaxies_tanh():