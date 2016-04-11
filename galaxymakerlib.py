# !/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
## 			This is a library of many functions to the code galaxymaker		##
## 						Created by Geferson Lucatelli 						##
##						lucatelli_asph_r2@virgophysics						##
##																			##
##		(original version 03.16.2016) (lastest version 04.06.2016)			##
##############################################################################

##############################################################################
##  IMPORTING SOME USEFULL MODULES
##############################################################################
from __future__ import division
import numpy as np
from numpy import meshgrid, arange, log, sqrt, tanh, tan, arctan, arctanh, arctan2
import matplotlib.pyplot as plt
import scipy as sp
import pylab as pl
import math
import pyfits as pf
import os
import requests
import sys


##############################################################################
##  SOME SURFACE BRIGHTNESS PROFILES
##############################################################################
'''
Some surface brightness and routines related to it, to create disk galaxies, elliptical galaxies and spiral galaxies.
'''
def grid(q,c,PA,gal_center=(0,0),N=4,M=4,size=255):
	'''
		Create an square meshgrid array (2D array), to create the galaixes' surface brightness.
		Retunrs an image array (pixels=size*size), an 1D arrays (x and y).
		Input Parameters:

		N,M 		spacial size of the image/arrays (for create a spiral structure, this values must lie approximately between 4 and 8)
		size		pixel scale size of the image (this is the width of the image)
		x,y 		1D arrays of the meshgrig
		gal_center 	the center of the galaxy
		r 			the radial grid from the center, or the radial profile of the image (used to create the image)
		c 			form of the generalizaed ellipse (c=0 is normal ellipse)
		delta		parameter to eliminate divisions by zero in the galaxy center
	'''
	delta=0.001
	x,y=meshgrid(arange(-M/2,M/2,M/size), arange(-N/2,N/2,N/size))
	x,y=rotation(PA,gal_center,x,y)
	r=((abs(x-gal_center[1])**(c+2.0)+((abs(y-gal_center[0]))/(q))**(c+2.0))**(1.0/(c+2.0))+delta)
	return r,x,y

def rscale(r_ef,N=4,size=255):
	'''
		A simple function to rescale the inputs parameters related to the efective radius of a given profile.
		returns the rescaled parameter (scalar)

		r_ef		the efective radius of any profile
					ex: rn, re, ...
	'''
	return r_ef*N/size

def rotation(PA, gal_center,x, y):
	#convert to radians
	t=PA*np.pi/180.
	return ((x-gal_center[1])*np.cos(t) + (y-gal_center[0])*np.sin(t), -(x-gal_center[1])*np.sin(t) + (y-gal_center[0])*np.cos(t))

def sersic_profile(In,rn,n,q,c):
	'''
		The Sersic Profile: used to create ellipticals and also spiral galaxies.
		The Sersic Profile is used in all functions bellow.
		returns an image array

		Input Paremeters:

		n 			Sersic index, control the concentration of galaxy light in relation to his center
					high n, high concentration of light
					low n, low concentration
					examples:
						n=1				exponential profile
						n=4				de vaucouleurs profile
						n=0.5 to n=3	disk/bar galaxies
						n=3 to 10		bulge/ellipticals galaxies in general
		In 			Sersic half light radius intensity
		rn 			radius containing the In [I(rn)=In]
		bn 			related quantitie with the sersic index
	'''
	r,x,y=grid(q,c,00)
	bn=1.9992*n-0.3271
	return In*np.exp(	-bn*(	(r/(rscale(rn))	)**(1.0/n) -1	)	)

def disk_profile(Id,rd,nd,q,c):
	'''
		Disk profile, to use with the spiral function
		returns an image array

		Id 			the disk half light radius intensity
		rd 			disk half light radius
		nd 			sersic index for the disk
	'''
	return sersic_profile(Id,rd,nd,q,c)

def bar_profile(Ib,rb,nb,q,cbar):
	'''
		Bar profile generates a square structure, something like a bar in the center. 
		This function is in test.
		In the bar profile, c can be more than 1
		It returns an image array


		Ib 			the bar half light radius intensity
		rb 			bar half light radius
		nb 			sersic index for the bar
	'''
	return sersic_profile(Ib,rb,nb,q,cbar)

def exponential_profile(Ie,re,q,c):
	'''
		Exponential profile, particular case of the sersic profile
		Returns an image array

		Ie 			the exponential half light radius intensity
		re 			exponential half light radius
		ne 			sersic index for the exponential: n=1
	'''
	ne=1.0
	return sersic_profile(Ie,re,ne,q,c)

def bulge_profile(Ibu,rbu,nbu,q,c):
	'''
		Bulge profile, an elliptical structure like, in the center of the image
		Returns an image array

		Ibu 		the bulge half light radius intensity
		rbu 		bulge half light radius
		nbu 		sersic index for the bulge
	'''
	return sersic_profile(Ibu,rbu,nbu,q,c)

def tan_spiral_profile(k,p,AA,NN,phi0,q):
	'''
		This function, generates a tangent spiral structure
		Returns an image array 

		k			the number of spiral arms
		p 			coefficient of proportionately of the spiral structure (the spiral is inversally proportionately to it parameter)
		NN 			the winding number, control the winding/tightening of the spiral arms
		phi0		control how may turns the spirals do
		BB 			related quantitie with phi0
		x,y 		are the arrays of the meshgrid (function grid()): is required because the use of arctan2(x,y)
		delta_tan	parameter to eliminate invalid values from the tangent function
		phi_r_tan	an angular function as function of r, used to create the spiral structure (spiral)
	'''
	c=0.0
	r,x,y=grid(q,c,30)
	delta_tan=0.01
	BB=(1.0)/(np.tanh(phi0/(2.0*NN)))
	phi_r_tan=2.0*NN*arctan(np.exp(AA/(r**1.0+delta_tan) )/BB)
	# I_exp=exponential_profile(1.0,20.0,q,c)
	I_exp=np.exp(-r/rscale(30.0))
	sp=p*np.cos(k*(arctan2(x,y)+phi_r_tan))
	spiral=((1.+sp)/(2))*I_exp
	# plot_two_profiles(spiral,I_exp)
	return spiral

def tanh_spiral_profile(k,p,AA,NN,phi0,q):
	'''
		This function, generates a hyperbolic tangent spiral structure
		Returns an image array

		k			the number of spiral arms
		p 			coefficient of proportionately of the spiral structure (the spiral is inversally proportionately to it parameter)
		NN 			the winding number, control the winding/tightening of the spiral arms
		phi0		control how may turns the spirals do
		BB 			related quantitie with phi0
		x,y 		are the arrays of the meshgrid (function grid()): is required because the use of arctan2(x,y)
		delta_tanh	parameter to erase invalid values from the hyperbolic tangent function
		phi_r_tanh	an angular function as function of r, used to create the spiral structure (spiral)
	'''
	c=0.0
	r,x,y=grid(q,c,30)
	delta_tanh=1.0
	BB=(1.0)/(tanh(phi0/(2.0*NN)))
	phi_r_tanh=2.0*NN*arctanh(np.exp(AA/(r**1.0+delta_tanh) )/BB)
	I_exp=np.exp(-r/rscale(30.0))
	sp=p*np.cos(k*(arctan2(x,y)+phi_r_tanh))
	spiral=((1.+sp)/(2))*I_exp
	# plot_two_profiles(spiral,I_exp)
	return spiral


##############################################################################
##  CREATE SYNTHETIC GALAXIES IMAGES: ELLIPTICALS AND SPIRALS
##############################################################################
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
'''
def spiral_galaxie(nsc,In,rn,n,ec,Ie,re,sct,k,p,AA,NN,phi0,barc,Ib,rb,nb,qbar,cbar,q,c):
	return nsc*sersic_profile(In,rn,n,q,c)+ec*exponential_profile(Ie,re,q,c)+sct*tan_spiral_profile(k,p,AA,NN,phi0,q)+barc*bar_profile(Ib,rb,nb,qbar,cbar)


##############################################################################
## DISPLAYING THE RESULTS
##############################################################################
def plot_save_image(image,number_name,nsc,In,rn,n,ec,Ie,re,sct,k,p,AA,NN,phi0,barc,Ib,rb,nb,q,c):
	fig=plt.figure()
	ax1=fig.add_subplot(2,2,1)
	ax1.imshow(((image))**0.2)
	ax2=fig.add_subplot(2,1,2)
	plt.xlabel('$r$ [pixels]')
	plt.ylabel('$I(r)$')
	ax2.plot(image[len(image)/2,len(image)/2:])
	ax3=fig.add_subplot(2,2,2)
	plt.axis('off')
	ax3.text(0.0, 1.00, r'Galaxy ID:$'+str(1+number_name)+'$', fontsize=12)
	#sersic/bulge parameters
	ax3.text(0.0, 0.92, r'$q$='+'$'+str(format(q,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.84, r'$r_n$='+'$'+str(format(rn,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.76, r'$I_n$='+'$'+str(format(In,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.68, r'$n$='+'$'+str(format(n,'.2f'))+'$', fontsize=12)
	#bar parameters
	ax3.text(0.0, 0.60, r'$r_b$='+'$'+str(format(rb,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.52, r'$I_b$='+'$'+str(format(Ib,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.44, r'$n_b$='+'$'+str(format(nb,'.2f'))+'$', fontsize=12)
	#exponential parameters
	ax3.text(0.0, 0.36, r'$r_e$='+'$'+str(format(re,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.28, r'$I_e$='+'$'+str(format(Ie,'.2f'))+'$', fontsize=12)
	#spiral parameters
	ax3.text(0.0, 0.20, r'$N$='+'$'+str(format(NN,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.12, r'$\phi_0$='+'$'+str(format(phi0,'.2f'))+'$', fontsize=12)
	ax3.text(-0.2, 0.04, r'$I(r)='+str(nsc)+'I_{ser}(r)+'+str(ec)+'I_{exp}(r)+'+str(barc)+'I_{bar}(r)+'+str(sct)+'\Theta(\phi(r))$',fontsize=12)
	plt.savefig('gal_'+str(1+number_name)+'.jpg',dpi=200)
	plt.show()
	# plt.clf()
	return

def plot_image(image):
	fig=plt.figure()
	ax1=fig.add_subplot(1,1,1)
	ax1.imshow(((image)))
	plt.show()
	return

def plot_two_profiles(profile1,profile2):
	'''
		Creates a plot comparing two profiles. It is useful to study how an exponential profile changes a spiral one.
		Retuns a graphic
	'''
	fig=plt.figure()
	ax=fig.add_subplot(1,1,1,axisbg='white')
	ax.plot(profile1[len(profile1)/2,len(profile1)/2:],'--k',color='red',label='$I_1$')
	ax.plot(profile2[len(profile2)/2,len(profile2)/2:],'--k',color='blue',label='$I_2$')
	ax.plot(profile1[len(profile1)/2,len(profile1)/2:]*profile2[len(profile2)/2,len(profile2)/2:],'--k',color='green',label='$I_1\\times I_2$')
	plt.legend(loc='uper left')
	plt.legend(loc=(0.80,0.70))
	plt.grid(True)
	plt.xlabel('$r$ [px]')
	plt.ylabel('$I(r)$')
	# plt.show()
	plt.savefig('gal_'+str(1+number_name)+'.jpg',dpi=200)
	plt.clf()
	return

def save_fits(image,number_name):
	'''
		Save to .fits the synthetic galaxie
		Returns a .fits file on the folder. 
	'''
	pf.writeto('gal_'+str(1+number_name)+'.fits', image,  clobber=1)
	print 'Gal', number_name+1, ' done'
	return