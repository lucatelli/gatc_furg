# !/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
## This is a main script that runs routines to create synthetic galaxies    ##
## 			Created by Geferson Lucatelli 			    ##
##			lucatelli_asph_r2@virgophysics 	 		    ##
##					  			            ##
##		(original version 03.16.2016) (lastest version 04.10.2016)  ##
##############################################################################

##############################################################################
##IMPORT SCRIPTS
##############################################################################
import galaxymakerlib
from galaxymakerlib import*
reload (galaxymakerlib)
import os
import requests
import sys

'''
Create synthetic galaxies, including in ellipticals and spirals components like bulges, disks and bars.
'''
'''
	The function create_galaxies() calls a specified profile, to create a synthetic galaxie. 
	returns an image array

		Some components that can be used to control the intensity of the profile in the galaxy are:
		
		ec	exponential contribution to the final galaxy intensity output
		
		buc	bulge contribution to the final galaxy intensity output
		
		dc 	disk contribution to the final galaxy intensity output
		
		barc	bar contribution to the final galaxy intensity output
		
		nsc	general sersic profile contribution to the final galaxy intensity output
		
		sct	tangent spiral contribution to the final galaxy intensity output
		
		scth	hyperbolic tangent spiral contribution to the final galaxy intensity output

		gal_center_ass 		include a profile in any other center. 
		It is useful to create assimetric galaxies.



		It can be created a simple synthetic galaxie in a more easy way simply calling a function.
'''
def create_galaxies():
	'''
		This function join some functions of galaxymakerlib to create a more complex sythetic galaxie, 
		including spiral arms with the use of normal tangent function. 
		See paper http://arxiv.org/abs/0908.0892
		Returns an array image.
	'''
	nber__galaxies=1
	'''
		A 'good' set of parameters to create galaxies with the normal and the hyperbolic tangent.
		There is one pint that must be observed: to create galaxies with the normal tangen, phi0 
		lies in [0.1,3.0], and NN in [1.0,4.0]. 
		However, with the hyperbolic tangent, these parameters lies on, for example phi0=[1.0,3.5], 
		NN=[>2*max(phi0),<4*max(phi0)]
	'''
	rn=[40.0]
	n=[5.0]
	# Fn=[10000]
	re=[10.0]
	Fn=np.linspace(8000.0,12000.0,nber__galaxies)
	sct=intensity_total_flux(1.0,30.0,Fn/20.)
	Ie=intensity_total_flux(1.0,re[0],Fn/500.)
	In=intensity_total_flux(n[0],rn[0],Fn)
	# rn=np.linspace(15.0,60.0,nber__galaxies)
	# n=np.linspace(2.0,6.0,nber__galaxies)
	# Ie=np.linspace(10.0,40.0,nber__galaxies)
	# re=np.linspace(15.0,35.,nber__galaxies)
	# sct=np.linspace(100.0,350.0,nber__galaxies)
	# NN=np.linspace(1.0,3.0,nber__galaxies)
	# phi0=np.linspace(0.1,2.4,nber__galaxies)
	# rb=np.linspace(5.0,15.0,nber__galaxies)
	# Ib=np.linspace(10.0,20.0,nber__galaxies)
	# q=np.linspace(0.5,1.0,nber__galaxies)
	# PA=np.linspace(0.0,360.0,nber__galaxies)
	# NN=np.linspace(8.0,16.0,nber__galaxies)
	# phi0=np.linspace(0.5,3.6,nber__galaxies)
	# In=[intensity_total_flux(n[0],rn[0],Fn[0])]
	NN=[2.0]
	phi0=[1.0]
	Ib=[15.0]
	rb=[10.0]
	q=[0.80]
	PA=[10.0]
	# #####################################################################################
	nsc=1.0
	ec=1.0
	k=2.0
	p=1.0
	AA=1.5
	barc=0.0
	nb=0.5
	qbar=0.2
	cbar=2.2
	c=0.0
	i=0
	gal_center=(0,0)
	gal_center_ass=(0,0)
	for j in range(nber__galaxies):
		values=(In[j],n[i],rn[i],Ie[i],re[i],sct[i],NN[i],phi0[i],rb[i],Ib[i],q[i],PA[i])
		gal=spiral_galaxie(nsc,In[j],rn[i],n[i],ec,Ie[i],re[i],sct[i],k,p,AA,NN[i],phi0[i],\
		barc,Ib[i],rb[i],nb,qbar,cbar,q[i],c,PA[i],gal_center,gal_center_ass)
		plot_save_image(gal,j,nsc,In[j],rn[i],n[i],ec,Ie[i],re[i],sct[i],k,p,AA,NN[i],phi0[i],\
		barc,Ib[i],rb[i],nb,qbar,cbar,q[i],c,PA[i])
		# w_values('gal',values)
		w_value('In_var',In[j])
		print 'gal_',j+1, ' done.'
		save_fits(gal,j)
		# plt.clf()
create_galaxies()
nber__galaxies=1
def call_mfmtk():
	for i in range(1,nber__galaxies+1):
		os.system("python path gal_"+str(i)+".fits psf.fits")
		print 'Gal ', str(i),' done.'
# call_mfmtk()
