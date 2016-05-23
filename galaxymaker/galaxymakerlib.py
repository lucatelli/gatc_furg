# !/usr/bin/env python
# -*- coding: utf-8 -*-
#####################################################################
## 	This is a library of many functions to the code galaxymaker	   ##
##	 				 Created by Geferson Lucatelli 			   	   ##
##					 lucatelli_asph_r2@virgophysics 	 		   ##
##					  			            					   ##
##	  (original version 03.16.2016) (lastest version 05.23.2016)   ##
#####################################################################

##############################################################################
##  IMPORTING SOME USEFULL MODULES
##############################################################################
from __future__ import division
import numpy as np
from numpy import meshgrid, arange, log, log10 ,sqrt, tanh, tan
from numpy import arctan, arctanh, arctan2,arcsinh
import matplotlib.pyplot as plt
import scipy as sp
from scipy.special import gammainc, gamma
import scipy.signal as sig
import pylab as pl
from pylab import where
import math
import pyfits as pf
import os
import requests
import sys
import scipy.ndimage as nd
import Image
import glob
from mfmtkutils import *
from morfometrykalibG2 import savitzky_golay_2d




##############################################################################
##  GENERAL TOOLS
##############################################################################
def gaussian2D(x0,y0,sigma, (M,N)):
	"""
	Return a 2D gaussian function.
	"""
	x,y = np.meshgrid(np.arange(N)-x0, np.arange(M)-y0)
	r2 = (x)**2 + (y)**2
	mask = 1./np.sqrt(2.*np.pi*sigma**2.) * np.exp(-r2/(2.*sigma**2.))
	return mask
def convolve(data,psf_name):
	return sig.fftconvolve(data, psf_name, 'same')

def rscale(r_ef,N=4.0,size=255):
	"""A simple function to rescale the inputs parameters related to 
	   the efective radius of a given profile.
       returns the rescaled parameter

       r_ef	the efective radius of any profile
					ex: rn, re, ...
	"""
	return r_ef*N/size

def rotation(PA,gal_center,x, y):
	"""Rotate an input image array. It can be used to modified 
	the position angle (PA).
	"""
	#convert to radians
	t=PA*np.pi/180.
	return ((x-rscale(gal_center[1]))*np.cos(t) + \
		(y-rscale(gal_center[0]))*np.sin(t),\
		-(x-rscale(gal_center[1]))*np.sin(t) + \
		(y-rscale(gal_center[0])*np.cos(t)))
def bn(n):
	"""Return the relation between the bn parameter with it's Sersic\
	index n. See Ciotti, L. and Bertin, G.; 1999. 
	"""
	return 1.9992*n-0.3271

def C1f(n):
	"""Returns the theoretical expression found for concentration C1 index 
	by Ferrari, F., Morfometryka, 2015. 
	"""
	return 2.91*((n/(32.44))**0.48)

def C2f(n):
	"""Returns the theoretical expression found for concentration C2 index 
	by Ferrari, F., Morfometryka, 2015. 
	"""
	return 1.39*((n/(13.36))**0.52)

def RPe(n,rn):
	"""Returns the theoretical expression found for petrosian radius, 
	by Ferrari, F., Morfometryka, 2015. 
	"""
	Rpm=5.8
	no=-1.11
	a=2.04
	al=0.8
	return rn*Rpm*((n-no)/a)*np.exp(-((n-no)/a)**al)

def Crn(n,aa=1.0/3.0):
	"""Returns the theoretical expression for the 'old' concentration index.
	See Graham, A.W. and Driver, S.P., 2005. 
	"""
	return (gammainc(2*n,bn(n)*((aa)**(1./n))))/(gammainc(2*n,bn(n)))
def mag(In,rn,n,r):
	"""Return the magnitude scale of a given light intensity profile.
	"""
	return +2.5*np.log10(In)+\
	((2.5*bn(n))/(np.log(10.0)))*((r[len(r)/2+1,len(r)/2:]/rn)**(1./n)-1)

# def eta_R(rn,n):
# 	R=RP(n,rn)
# 	x=bn(n)*((R/rn)**(1./n))
# 	return ((2*n)/(		(x**(2.*n))*np.exp(-x)		))*(gammainc(2*n,x))

##############################################################################
##  SOME SURFACE BRIGHTNESS PROFILES
##############################################################################
'''
Some surface brightness and routines related to it, to create disk galaxies, 
elliptical galaxies and spiral galaxies.
'''
def grid(q,c,PA,gal_center,N=4.0,M=4.0,size=255.0):
	"""Create an square meshgrid array (2D array), to create the galaixes' 
	surface brightness.
	Retunrs an image array (pixels=size*size), an 1D arrays (x and y).
    
    Input Parameters:

		N,M 		spacial size of the image/arrays (for create a spiral 
		structure, this values must lie approximately between 4 and 8)
		
		size		pixel scale size of the image (this is the width of the image)
		
		x,y 		1D arrays of the meshgrig
		
		gal_center 	the center of the galaxy
		
		r 		the radial grid from the center, or the radial profile 
		of the image (used to create the image)
		
		c 		form of the generalizaed ellipse (c=0 is normal ellipse)
		
		delta		parameter to eliminate divisions by zero in the galaxy center
	"""
	delta=0.001
	x,y=meshgrid(arange(-M/2,M/2,M/size), arange(-N/2,N/2,N/size))
	x,y=rotation(PA,gal_center,x,y)
	r=((abs(x-rscale(gal_center[1]))**(c+2.0)+\
		((abs(y-rscale(gal_center[0])))/(q))**(c+2.0))**(1.0/(c+2.0))+delta)
	return r,x,y

def radial_profile(data, center,q,c):
	"""Return the radial light profile (isophotes of light) of a given 
	surface brightness.
	"""
	y, x = np.indices((data.shape))
	# R = np.sqrt((x - center[0])**2 + (y - center[1])**2)
	R=(abs(x-center[1])**(c+2.0)+((abs(y-center[0]))/(q))**(c+2.0))**(1.0/(c+2.0))
	R = R.astype(np.int)
	tbin = np.bincount(R.ravel(), data.ravel())
	nr = np.bincount(R.ravel())
	radialprofile = tbin / 5.
	return radialprofile

def sersic_profile(In,rn,n,q,c,gal_center):
	"""The Sersic Profile: used to create ellipticals and also 
	spiral galaxies.
	The Sersic Profile is used in all functions bellow.
	returns an image array

	Input Paremeters:

		n 	Sersic index, control the concentration of 
		galaxy light in relation to his center
		
				high n, high concentration of light
				low n, low concentration
				examples:
				n=1		exponential profile
				n=4		de vaucouleurs profile
				n=0.5 to n=3	disk/bar galaxies
				n=3 to 10	bulge/ellipticals galaxies in general
				
		In 	Sersic half light radius intensity
		
		rn 	radius containing the In [I(rn)=In]
		
		bn 	related quantitie with the sersic index
	"""
	r,x,y=grid(q,c,0.0,gal_center)
	bn=1.9992*n-0.3271
	# F=(In*rn**2.0*2*np.pi*n*np.exp(bn))/(bn**(2*n))*gammainc(2*n,	(0.8*r.max()/rn)**(1.0/n)	)
	# print 'In=', In, ' F=',F
	return In*np.exp(	-bn*(	(r/(rscale(rn))	)**(1.0/n) -1	))
def luminosity_R(R_rat,In,n):
	return (In*(R_rat)**2.0*2*np.pi*n*np.exp(bn(n) 	))/(bn(n)**(2*n))*gamma(2*n)
def bulge_disk_component(n,rn,In,re,BT,alpha=1.0,beta=1.0):
	'''
		Convert the given flux variable quantitie, to the half light intensity In=(rn)
		Return the number In. 
	'''
	bn=1.9992*n-0.3271
	# In=((bn**2.0)*np.exp(-bn)/((rn)**2.0*2*np.pi*n))*Fn*(1./gamma(2*n))
	L_BT=(In*rn**2.0*2*np.pi*n*np.exp(bn))/(bn**(2*n))*gamma(2*n)
	L_DT=(alpha/beta)*(L_BT/BT)*(1.-BT)
	L_T=L_BT+L_DT
	BD=L_BT/L_DT
	Ie=((L_DT)/(2*np.pi*re**2.0*gamma(2)*np.exp(1.9992-0.3271)))*((1.9992-0.3271)**2.0)
	# print Ie
	return Ie, L_BT,L_DT,L_T,BD
def bulge_disk(BD,n,rn,re,LT,alpha=1.0,beta=1.0):
	L_BT=(BD*LT)/(beta+alpha*BD)
	L_DT=(LT-alpha*L_BT)/(beta)
	BT=1-(beta*L_BT)/(LT*BD)
	In=(bn(n)**(2.0*n)*L_BT)/(2.0*np.pi*n*((rn)**2.0)*np.exp(bn(n))*gamma(2.0*n))
	mu_n=-2.5*np.log10(In)
	mean_mu_n=mu_n-2.5*np.log10(	gamma(2.0*n)*(n*np.exp(bn(n))/(bn(n)**(2.0*n))	))
	Ie=((L_DT)/(2*np.pi*(re)**2.0*gamma(2)*np.exp(1.9992-0.3271)))*((1.9992-0.3271)**2.0)
	mu_e=-2.5*np.log10(Ie)
	mean_mu_e=mu_e-\
	2.5*np.log10(	gamma(2.0*1.0)*(1.0*np.exp(bn(1.0))/(bn(1.0)**(2.0*1.0))	))
	# print 'LBT=',L_BT,' LDT=',L_DT,'  BT=', BT,'  In=',In,'  Ie=', Ie
	return L_BT,L_DT,BT,In,Ie,mu_n,mean_mu_n,mu_e,mean_mu_e
def Inn(n,rn,L_BT):
	return (bn(n)**(2.0*n)*L_BT)/(2.0*np.pi*n*(rn**2.0)*np.exp(bn(n))*gamma(2.0*n))
# def disk_component(L_BT,BT,alpha=1.0,beta=1.0):
# 	'''
# 		Return the disk-to-total L_DT component of the galaxy luminoity, of a 
# 		given bulge-to-total component

# 		alpha	coefficient of proportionality of L_BT
# 		beta 	coefficient of proportionality of L_DT
# 		BT 		bulge-to-total luminosity
# 		L_BT	total bulge luminosity

# 	'''
# 	L_DT=(alpha/beta)*(L_BT/BT)*(1.-BT)
# 	Ie=((L_DT)/(2*pi*re**2.0*gamma(2)*np.exp(1.9992-0.3271)))*((1.9992-0.3271)**2.0)
# 	return Ie

def disk_profile(Id,rd,nd,q,c,gal_center):
	'''
		Disk profile, to use with the spiral function
		returns an image array

		Id 	the disk half light radius intensity
		
		rd 	disk half light radius
		
		nd 	sersic index for the disk
	'''
	return sersic_profile(Id,rd,nd,q,c,gal_center)

def bar_profile(Ib,rb,nb,q,cbar,gal_center):
	'''
		Bar profile generates a square structure, something like a 
		bar in the center. 
		
		This function is in test. In the bar profile, c can be more than 1
		Returns an image array


		Ib 	the bar half light radius intensity
		
		rb 	bar half light radius
		
		nb 	sersic index for the bar
	'''
	return sersic_profile(Ib,rb,nb,q,cbar,gal_center)

def exponential_profile(Ie,re,q,c,gal_center):
	'''
		Exponential profile, particular case of the sersic profile
		Returns an image array

		Ie 	the exponential half light radius intensity
		
		re 	exponential half light radius
		
		ne 	sersic index for the exponential: n=1
	'''
	ne=1.0
	return sersic_profile(Ie,re,ne,q,c,gal_center)

def bulge_profile(Ibu,rbu,nbu,q,c,gal_center):
	'''
		Bulge profile, an elliptical structure like, in the center of the image
		Returns an image array

		Ibu 	the bulge half light radius intensity
		rbu 	bulge half light radius
		nbu 	sersic index for the bulge
	'''
	return sersic_profile(Ibu,rbu,nbu,q,c,gal_center)

def tan_spiral_profile(k,spo,AA,NN,phi0,q,PA,gal_center,re,Ie):
	'''
		This function, generates a tangent spiral structure
		Returns an image array 

		k	the number of spiral arms
		
		p 	coefficient of proportionately of the spiral structure 
		(the spiral is inversally proportionately to it parameter)
		
		NN 	the winding number, control the winding/tightening of 
		the spiral arms
		
		phi0	control how may turns the spirals do
		BB 	related quantitie with phi0
		
		x,y 	are the arrays of the meshgrid (function grid()): 
		is required because the use of arctan2(x,y)
		
		delta_tan	parameter to eliminate invalid values from 
		the tangent function
		
		phi_r_tan	an angular function as function of r, used to 
		create the spiral structure (spiral)
	'''
	c=0.0
	r,x,y=grid(q,c,PA,gal_center)
	delta_tan=0.01
	BB=(1.0)/(np.tanh(phi0/(2.0*NN)))
	phi_r_tan=2.0*NN*arctan(np.exp(AA/(r**1.0+delta_tan) )/BB)
	# I_exp=exponential_profile(1.0,20.0,q,c)
	I_exp=np.exp(-r/rscale(re*1.0))
	sp=np.cos(k*(arctan2(x,y)+phi_r_tan))
	spiral=spo*Ie*((1.+sp)/(2))*I_exp
	# plot_two_profiles(spiral,I_exp)
	return spiral

def tanh_spiral_profile(k,p,AA,NN,phi0,q,PA,gal_center,re,Ie):
	'''
		This function, generates a hyperbolic tangent spiral structure
		Returns an image array

		k	the number of spiral arms
		
		p 	coefficient of proportionately of the spiral structure 
		(the spiral is inversally proportionately to it parameter)
		
		NN 	the winding number, control the winding/tightening of 
		the spiral arms
		
		phi0	control how may turns the spirals do
		
		BB 	related quantitie with phi0
		
		x,y 	are the arrays of the meshgrid (function grid()): 
		is required because the use of arctan2(x,y)
		
		delta_tanh	parameter to erase invalid values from 
		the hyperbolic tangent function
		
		phi_r_tanh	an angular function as function of r, 
		used to create the spiral structure (spiral)
	'''
	c=0.0
	r,x,y=grid(q,c,PA,gal_center)
	delta_tanh=1.0
	BB=(1.0)/(tanh(phi0/(2.0*NN)))
	phi_r_tanh=2.0*NN*arctanh(np.exp(AA/(r**1.0+delta_tanh) )/BB)
	I_exp=np.exp(-r/rscale(re))
	sp=p*np.cos(k*(arctan2(x,y)+phi_r_tanh))
	spiral=Ie*((1.+sp)/(2))*I_exp
	# plot_two_profiles(spiral,I_exp)
	return spiral


##############################################################################
##  CREATE SYNTHETIC GALAXIES IMAGES: ELLIPTICALS AND SPIRALS
##############################################################################
'''
Create synthetic galaxies, including in ellipticals and spirals components 
like bulges, disks and bars.
'''
'''
	The function galaxie() calls a specified profile, to create a synthetic \
	galaxie. 
	
	returns an image array

		Some components that can be used to control the intensity of 
		the profile in the galaxy are:
		
		ec	exponential contribution to the final galaxy intensity output
		
		buc	bulge contribution to the final galaxy intensity output
		
		dc 	disk contribution to the final galaxy intensity output
		
		barc	bar contribution to the final galaxy intensity output
		
		nsc	general sersic profile contribution to the final galaxy intensity \
		output
		
		sct	tangent spiral contribution to the final galaxy intensity output
		
		scth	hyperbolic tangent spiral contribution to the final galaxy \
		intensity output
'''

def bar(barc,Ib,rb,nb,qbar,cbar,gal_center):
	return barc*bar_profile(Ib,rb,nb,qbar,cbar,gal_center)

def ell_gal(n,In,rn,Ie,re,ec,nsc,qb,qd,c,PA,gal_center):
	Iser=sersic_profile(In,rn,n,qb,c,gal_center)
	Iexp=exponential_profile(Ie,re,qd,c,gal_center)
	return nsc*Iser+ec*Iexp, Iser,Iexp

##############################################################################
## DISPLAYING THE RESULTS
##############################################################################
def plot_save_gal(image,sersic,disk,spiral,n,nn,In,rn,Ie,re,qb,qd,c,PA,BT,BD,\
	L_BT,L_DT,L_T,nsc,ec,k,spo,AA,NN,phi0,number_name,path,SPIRAL,DISK):
	
	'''
	plot_save_gal() gets the input parameters and plot all the related \
	calculation in a figure, composed by:
		1 - the model galaxy output;
		2 - table with the input and output parameters;
		3 - a line plot of the surface brightness;
		4 - a galaxy polar plot. Code performed by Ferrari, F. on Morfometryka (2015). \
		See polarim() function.

	The variable SPIRAL is to make easy the plots when there is or not \
	a spiral component. If it is True, plots related to the spiral will \
	be computed, otherwise, it will not.
	
	Save a .jpg (or a .png or .pdf) file.
	'''

	fig=plt.figure()
	fig.set_size_inches(15.0, 10.0)
	ax1=fig.add_subplot(2,2,1)
	plt.title('Galaxia Sintetica', fontsize=10)
	ax1.imshow((np.arcsinh(image)))


	ax2=fig.add_subplot(2,2,3)
	plt.title('Perfis de Brilho', fontsize=10)
	plt.xlabel('$r$ [pixels]')
	plt.ylabel('$\mu(r)=2.5 \log_{	10}[I(r)]$')
	plt.xlim(1,100)
	plt.ylim(-1,10)
	ax2.plot(2.5*np.log10(sersic[len(sersic)/2,len(sersic)/2:]),'o',\
		color='red',label='$\mu_{bojo}$')
	if SPIRAL is True:
		plt.ylim((2.5*np.log10(spiral[len(spiral)/2,len(spiral)/2:])).min(),\
			(2.5*np.log10(image[len(image)/2+1,len(image)/2:])).max())
		ax2.plot(2.5*np.log10(spiral[len(spiral)/2,len(spiral)/2:]),'o',\
			color='green',label='$\mu_{spiral}$')
		
		if DISK is True:
			ax2.plot(2.5*np.log10(disk[len(disk)/2,len(disk)/2:]),'o',\
				color='blue',label='$\mu_{disco}$')
			ax2.plot(2.5*np.log10(image[len(image)/2,len(image)/2:]),'o',\
				color='black',label='$\mu_{bojo}+\mu_{disco}+\mu_{spiral}$')
		
		else:
			ax2.plot(2.5*np.log10(image[len(image)/2,len(image)/2:]),'o',\
				color='black',label='$\mu_{bojo}+\mu_{spiral}$')
	
	if DISK is True:
		ax2.plot(2.5*np.log10(disk[len(disk)/2,len(disk)/2:]),'o',\
			color='blue',label='$\mu_{disco}$')
		
		if SPIRAL is True:
			ax2.plot(2.5*np.log10(spiral[len(spiral)/2,len(spiral)/2:]),'o',\
				color='green',label='$\mu_{spiral}$')
			ax2.plot(2.5*np.log10(image[len(image)/2,len(image)/2:]),'o',\
				color='black',label='$\mu_{bojo}+\mu_{	disco}+\mu_{spiral}$')
		
		else:
			ax2.plot(2.5*np.log10(image[len(image)/2,len(image)/2:]),'o',\
				color='black',label='$\mu_{bojo}+\mu_{disco}$')
	
	plt.legend(loc='upper right',fancybox=True, shadow=True)
	plt.grid()


	ax4=fig.add_subplot(2,2,4)
	plt.title('Polar Plot', fontsize=10)
	ax4 = plt.gca()
	ax4.imshow(np.log((polarim(image))[:150]))
	ax4.invert_yaxis()
	# ax2.text(1.0*rn, plt.ylim()[0] , r'$\downarrow$', color='red', fontsize=25)


	ax3=fig.add_subplot(2,2,2)
	plt.axis('off')
	ax3.text(0.0, 1.10, r'Galaxia ID:$'+str(number_name)+'$', fontsize=12)
	#sersic/bulge/disk parameters
	
	if nn[1] is 'n1D':
		ax3.text(0.0, 1.00, r'$n_{	1D}$='+'$'+str(format(n,'.2f'))+'$',\
		 fontsize=12)
	
	if nn[1] is 'n2D':
		ax3.text(0.0, 1.00, r'$n_{	2D}$='+'$'+str(format(n,'.2f'))+'$',\
		 fontsize=12)
	
	if nn[1] is 'n':
		ax3.text(0.0, 1.00, r'$n$='+'$'+str(format(n,'.2f'))+'$',\
		 fontsize=12)
	# ax3.text(0.0, 1.00, r'$n_{	1D}$='+'$'+str(format(n,'.2f'))+'$', fontsize=12)
	
	ax3.text(0.0, 0.92, r'$I_n$='+'$'+str(format(In,'.2f'))+'$', \
		fontsize=12)
	ax3.text(0.0, 0.84, r'$R_n$='+'$'+str(format(rn,'.2f'))+'$', \
		fontsize=12)
	
	if DISK is True:
		ax3.text(0.0, 0.76, r'$I_d$='+'$'+str(format(Ie,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.0, 0.68, r'$R_d$='+'$'+str(format(re,'.2f'))+'$', \
			fontsize=12)
	else:
		ax3.text(0.0, 0.76, r'$I_d$='+'$'+str(format(0.00,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.0, 0.68, r'$R_d$='+'$'+str(format(0.00,'.2f'))+'$', \
			fontsize=12)

	#general ellipsis parameters
	
	if DISK is True:
		ax3.text(0.0, 0.60, r'$q_b$='+'$'+str(format(qb,'.2f'))+'$'+',  \
			$q_d$='+'$'+str(format(qd,'.2f'))+'$', fontsize=12)

	else:
		ax3.text(0.0, 0.60, r'$q_b$='+'$'+str(format(qb,'.2f'))+'$'+',  \
			$q_d$='+'$'+str(format(0.00,'.2f'))+'$', fontsize=12)
	
	ax3.text(0.0, 0.52, r'$c$='+'$'+str(format(c,'.2f'))+'$', fontsize=12)
	ax3.text(0.0, 0.44, r'$PA$='+'$'+str(format(PA,'.2f'))+'$', fontsize=12)
	#Output Gal Quantities
	ax3.text(0.0, 0.36, r'$\xi_{BT}$='+'$'+str(format(BT,'.2f'))+'$', \
		fontsize=12)
	
	if DISK is True:
		ax3.text(0.0, 0.28, r'$\xi_{BD}$='+'$'+str(format(BD,'.2f'))+'$', \
			fontsize=12)
	
	else:
		ax3.text(0.0, 0.28, r'$\xi_{BD}$='+'$'+'nan'+'$', fontsize=12)
	
	ax3.text(0.0, 0.20, r'$\mathscr{L}_{BT}$='+'$'+str(format(L_BT,'.2f'))+\
		'$', fontsize=12)
	
	if DISK is True:
		ax3.text(0.0, 0.12, r'$\mathscr{L}_{DT}$='+'$'+\
			str(format(L_DT,'.2f'))+'$', fontsize=12)
	
	else:
		ax3.text(0.0, 0.12, r'$\mathscr{L}_{DT}$='+'$'+\
			str(format(0.00,'.2f'))+'$', fontsize=12)
	
	ax3.text(0.0, 0.04, r'$\mathscr{L}_{T}$='+'$'+str(format(L_T,'.2f'))+'$',\
	 fontsize=12)
	
	if SPIRAL is True:
		ax3.text(0.5, 1.00, r'$N$='+'$'+str(format(NN,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.5, 0.92, r'$\phi_0$='+'$'+str(format(phi0,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.5, 0.84, r'$A$='+'$'+str(format(AA,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.5, 0.76, r'$\kappa$='+'$'+str(format(k,'.2f'))+'$', \
			fontsize=12)
		ax3.text(0.5, 0.68, r'$s_{p0}$='+'$'+str(format(spo,'.2f'))+'$', \
			fontsize=12)
	# ax3.text(-0.2, 0.04, r'$I(r)='+str(format(nsc,'.2f'))+'I_{bul}(r)+'+\
		# str(ec)+'I_{disk}(r)$',fontsize=12)
	'''Save the results'''
	plt.savefig(path+'jpg/gal_'+str(number_name)+'.jpg',dpi=200)
	# plt.show()
	plt.clf()
	
def plot_grad(image,sersic,disk,spiral,number_name,path,SPIRAL):
	
	'''
	plot_grad() takes the light gradient field from the polar decomposition of \
	the galaxy, function polarim(). Code performed by Ferrari, F. on \
	Morfometryka (2015). See zavitzky_golay_2d() and grad(). 

	It save a image file with the variation of the mean gradient field as \
	a function of the galactocentric radius of the galactocentric. 
	'''

	fig=plt.figure()
	g_BD=grad(polarim(sersic+disk)[:120])
	axx1=fig.add_subplot(1,1,1)
	fig.set_size_inches(10.0, 6.0)
	axx1.plot((g_BD),'--o',color='blue',label='bojo+disco')
	if SPIRAL is True:
		g_sp_BD=grad(polarim(spiral+disk+sersic)[:120])
		axx1.plot(g_sp_BD,'--o',color='red',label='spiral+bojo+disco')
	plt.title('Galaxia ID:'+str(number_name), fontsize=10)
	plt.xlabel('$r$ [pixels]')
	plt.ylabel('Gradiente')
	plt.legend(loc='upper right',fancybox=True, shadow=True)
	plt.grid()
	plt.savefig(path+'jpg/grad/grad_'+str(number_name)+'.jpg',dpi=200)
	plt.clf()
def plot_image(image_array):
	'''
		Fast way to plot an image. 
		Return an imshow() array.
	'''
	fig=plt.figure()
	ax1=fig.add_subplot(1,1,1)
	ax1.imshow(((image_array))**0.2)
	plt.show()
	return

def plot_two_profiles(profile1,profile2):
	'''
		Creates a plot comparing two profiles. It is useful to study how an 
		exponential profile changes a spiral one.
		Retuns a graphic
	'''
	fig=plt.figure()
	ax=fig.add_subplot(1,1,1,axisbg='white')
	ax.plot(profile1[len(profile1)/2,len(profile1)/2:],'--k',\
		color='red',label='$I_1$')
	ax.plot(profile2[len(profile2)/2,len(profile2)/2:],'--k',\
		color='blue',label='$I_2$')
	ax.plot(profile1[len(profile1)/2,\
		len(profile1)/2:]*profile2[len(profile2)/2,len(profile2)/2:],'--k',\
		color='green',label='$I_1\\times I_2$')
	plt.legend(loc='uper left')
	plt.legend(loc=(0.80,0.70))
	plt.grid(True)
	plt.xlabel('$r$ [px]')
	plt.ylabel('$I(r)$')
	# plt.show()
	plt.savefig('gal_'+str(number_name)+'.jpg',dpi=200)
	plt.clf()
	return

def save_fits(image,number_name,name,path):
	'''
		Save to .fits the synthetic galaxie
		Returns a .fits file on the folder. 
	'''
	if name is None:
		if number_name<=9:
			pf.writeto(path+'fits/gal_0'+str(number_name)+'.fits', image,  \
				clobber=1)
		else:
			pf.writeto(path+'fits/gal_'+str(number_name)+'.fits', image,  \
				clobber=1)
	else:
		pf.writeto(path+'fits/'+name+'.fits', image,  clobber=1)
	# print 'Gal', number_name+1, ' done'
	return

def save_fits_here(image,number_name,name):
	'''
		Save to .fits the synthetic galaxie
		Returns a .fits file on the folder. 
	'''
	if name is None:
		if number_name<=9:
			pf.writeto('gal_0'+str(number_name)+'.fits', image,  clobber=1)
		else:
			pf.writeto('gal_'+str(number_name)+'.fits', image,  clobber=1)
	else:
		pf.writeto(name+'.fits', image,  clobber=1)
	# print 'Gal', number_name+1, ' done'
	return


def w_value(file_name,p_value,path):
	'''
		Write on a file (line format) the needed input parameter given, \
		of each galaxie. 

		Save a file.dat on the folder.
	'''
	with open(path+file_name+'_values.dat', 'a') as f:
			f.write(str(p_value))
			f.write('')

def w_values(file_name,p_value,path):
	'''
		Write on a file the input parameters given, for each galaxies. \
		Each galaxy it 
		is a columm in the file.

		Save a file.dat on the folder.
	'''
	for i in p_value:
		with open(path+file_name+'_values.dat', 'a') as f:
			f.write(str(i))
			f.write(',')
	with open(path+file_name+'_values.dat', 'a') as f:
		f.write('\n')


##############################################################################
## PARAMETERS EVOLUTION
##############################################################################

def plots(var1,var2,labels,plot_config,path):
	""" Plots a wanted pair of variables, var1 and var2.

		var1: It is an input variable: for example, a gml one,
		or a mfmtk one.

		var2: It is the other variable. 

		plot_config: It is a simple configuration for the plot. For example, 
		can be some values of the variables that are kept constant in the 
		plot. The form of plot_config is: 
		plot_config=[['constant variable 1 string name',variable value,'LaTeX \
		code for this variable'],...]
		
		Example: 
		plot_config=[['qb',qb[0],'$q_b$'],['qd',qd,'$q_d$'],['n',n[0],'$n$'],\
		['BD',BD[0],'$\\xi_{	BD}$']]

	"""
	fig = plt.gcf()
	plt.plot(var1,var2,'--o',color='red')
	if labels[1] is 'C_1':
		plt.axhline(y=C1f(1.0), xmin=0, xmax=1, hold=None,label='Disco Puro',\
			color='purple')
		plt.axhline(y=C1f(6.0), xmin=0, xmax=1, hold=None,label='Bojo Puro' ,\
			color='brown')
	if labels[1] is 'C_2':
		plt.axhline(y=C2f(1.0), xmin=0, xmax=1, hold=None,label='Disco Puro',\
			color='purple')
		plt.axhline(y=C2f(6.0), xmin=0, xmax=1, hold=None,label='Bojo Puro' ,\
			color='brown')
	if labels[0] is 'C_1':
		plt.axvline(x=C1f(1.0), ymin=0, ymax=1, hold=None,label='Disco Puro',\
			color='purple')
		plt.axvline(x=C1f(6.0), ymin=0, ymax=1, hold=None,label='Bojo Puro' ,\
			color='brown')
	if labels[0] is 'C_2':
		plt.axvline(x=C2f(1.0), ymin=0, ymax=1, hold=None,label='Disco Puro',\
			color='purple')
		plt.axvline(x=C2f(6.0), ymin=0, ymax=1, hold=None,label='Bojo Puro' ,\
			color='brown')
	plt.title('$'\
			  +str(plot_config[0][2])+'='+str(plot_config[0][1][0])+'$, $'\
			  +str(plot_config[1][2])+'='+str(plot_config[1][1][0])+'$, $'\
			  +str(plot_config[2][2])+'='+str(plot_config[2][1][0])+'$, $'\
			  +str(plot_config[3][2])+'='+str(plot_config[3][1][0])+'$, $'\
			  +str(plot_config[4][2])+'='+str(plot_config[4][1][0])+'$',\
			  x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.92,x=0.70 , fontsize=13)
	plt.xlabel('$'+str(labels[0][0:])+'$',fontsize=20)
	plt.ylabel('$'+str(labels[1][0:])+'$',fontsize=20)
	plt.grid()
	fig.set_size_inches(12.0, 7.0)
	plt.legend(loc='best',fancybox=True, shadow=True)
	plt.legend(loc=(0.89,0.65))
	plt.savefig(path+'plots/'+str(labels[0][0:])+'_'+labels[1]+'.pdf')
	plt.savefig(path+'plots/'+str(labels[0][0:])+'_'+labels[1]+'.jpg')
	# plt.show()
	print 'Plot', labels[0], ' versus ', labels[1], ' done!'
	plt.clf()


def compare_n_R(path,path2,what_n,mfmtk_cat1,glm_cat1,mfmtk_cat2,glm_cat2):
	'''
	Forma de path2 -->> path2=[ 'diretorio', ' qual n eh, modelos com 1D ou 2D',\
	'Latex code' ]
	'''
	''' Para cada galaxia, verificar a relacao de (n1D, n2D)  do bojo+disco 
		com as medidas mormometricas, comparando com os mesmo (n1D,n2D) 
		de somente bojo com as medidas morfometricas deste.

		n1D_eff1: indice de sersic 1D effetivo do modelo bojo+disco
		n2D_eff1: indice de sersic 2D effetivo do modelo bojo+disco

		n1D_eff2: indice de sersic 1D efetivo do modelo bojo criado com
				  o indice de sersic n1D_eff1, espera-se que ambos sejam
				  iguais.

		n2D_eff2: indice de Sersic 2D efetivo do modelo bojo criaco com
				  o indice de Sersic n2D_eff1, espera-se que ambos sejam
				  iguais.
	'''
	fig = plt.gcf()
	print '--------------- LOADING RADIUS RATIOS-------------------'
	''' Lendo os Raios de luz do modelo bojo+disco e os indices de sersic efetivos \
	n_ef -->> n_d+n_b  para este modelo
	'''
	R10=mfmtk_cat1.param_selection(['R10'],column_dict)
	R20=mfmtk_cat1.param_selection(['R20'],column_dict)
	R30=mfmtk_cat1.param_selection(['R30'],column_dict)
	R40=mfmtk_cat1.param_selection(['R40'],column_dict)
	R50=mfmtk_cat1.param_selection(['R50'],column_dict)
	R60=mfmtk_cat1.param_selection(['R60'],column_dict)
	R70=mfmtk_cat1.param_selection(['R70'],column_dict)
	R80=mfmtk_cat1.param_selection(['R80'],column_dict)
	R90=mfmtk_cat1.param_selection(['R90'],column_dict)
	n1D_eff_1=mfmtk_cat1.param_selection(['nFit1D'],column_dict)
	n2D_eff_1=mfmtk_cat1.param_selection(['nFit2D'],column_dict)
	In2D=mfmtk_cat1.param_selection(['InFit2D'],column_dict)
	nbd=glm_cat1.param_selection(['n'],column_dict_glm)
	''' Lendo os Raios de Luz do modelo bojo, com indice de sersic efetivo \
	n_ef (obtidos de n1D ou n2D)
	'''
	Ref10=mfmtk_cat2.param_selection(['R10'],column_dict)
	Ref20=mfmtk_cat2.param_selection(['R20'],column_dict)
	Ref30=mfmtk_cat2.param_selection(['R30'],column_dict)
	Ref40=mfmtk_cat2.param_selection(['R40'],column_dict)
	Ref50=mfmtk_cat2.param_selection(['R50'],column_dict)
	Ref60=mfmtk_cat2.param_selection(['R60'],column_dict)
	Ref70=mfmtk_cat2.param_selection(['R70'],column_dict)
	Ref80=mfmtk_cat2.param_selection(['R80'],column_dict)
	Ref90=mfmtk_cat2.param_selection(['R90'],column_dict)
	n1D_eff_2=mfmtk_cat2.param_selection(['nFit1D'],column_dict)
	n2D_eff_2=mfmtk_cat2.param_selection(['nFit2D'],column_dict)
	In2Def=mfmtk_cat2.param_selection(['InFit2D'],column_dict)
	nb=glm_cat2.param_selection(['n'],column_dict_glm)

	R21=np.asarray(R20)/np.asarray(R10)
	R32=np.asarray(R30)/np.asarray(R20)
	R43=np.asarray(R40)/np.asarray(R30)
	R54=np.asarray(R50)/np.asarray(R40)
	R65=np.asarray(R60)/np.asarray(R50)
	R76=np.asarray(R70)/np.asarray(R60)
	R87=np.asarray(R80)/np.asarray(R70)
	R98=np.asarray(R90)/np.asarray(R80)
	Ref21=np.asarray(Ref20)/np.asarray(Ref10)
	Ref32=np.asarray(Ref30)/np.asarray(Ref20)
	Ref43=np.asarray(Ref40)/np.asarray(Ref30)
	Ref54=np.asarray(Ref50)/np.asarray(Ref40)
	Ref65=np.asarray(Ref60)/np.asarray(Ref50)
	Ref76=np.asarray(Ref70)/np.asarray(Ref60)
	Ref87=np.asarray(Ref80)/np.asarray(Ref70)
	Ref98=np.asarray(Ref90)/np.asarray(Ref80)
	'''Call a function to obtain the luminosities in these radius.'''
	print '----- CALCULATING THE LUMINOSITIES AT THE MANY RADIUS -----'
	L20=luminosity_R(np.asarray(R20),np.asarray(In2D),np.asarray(n2D_eff_1))
	L40=luminosity_R(np.asarray(R40),np.asarray(In2D),np.asarray(n2D_eff_1))
	L60=luminosity_R(np.asarray(R60),np.asarray(In2D),np.asarray(n2D_eff_1))
	L80=luminosity_R(np.asarray(R80),np.asarray(In2D),np.asarray(n2D_eff_1))
	LR20=L20/np.asarray(R20)**2.0
	LR40=L40/np.asarray(R40)**2.0
	LR60=L60/np.asarray(R60)**2.0
	LR80=L80/np.asarray(R80)**2.0
	print '--------------- PLOTING RESULTS -------------------'
	plt.plot(nbd,LR20,'--o',color='red',label='$L_{	20}/R_{	20}^2$')
	plt.plot(nbd,LR40,'--o',color='blue',label='$L_{	40}/R_{	40}^2$')
	plt.plot(nbd,LR60,'--o',color='green',label='$L_{	60}/R_{	60}^2$')
	plt.plot(nbd,LR80,'--o',color='black',label='$L_{	80}/R_{	80}^2$')
	# plt.xlabel('Indice de Sersic efetivo $n_{	ef}$: perfil bojo+disco.')
	plt.xlabel('$n$')
	plt.ylabel('Razao luminosidade raio: $L_i/R_i^2$')
	plt.grid()
	fig.set_size_inches(12.0, 7.0)
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.title('Modelos bojo+disco.',x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.93,x=0.70 , fontsize=13)
	plt.savefig(path2+'Li_Ri_n2D.jpg',dpi=300)
	plt.savefig(path2+'Li_Ri_n2D.pdf')
	# plt.show()
	plt.clf()
	print '--------------- 20% -------------------'
	Lef20=luminosity_R(np.asarray(Ref20),np.asarray(In2Def),\
		np.asarray(n2D_eff_2))
	Lef40=luminosity_R(np.asarray(Ref40),np.asarray(In2Def),\
		np.asarray(n2D_eff_2))
	Lef60=luminosity_R(np.asarray(Ref60),np.asarray(In2Def),\
		np.asarray(n2D_eff_2))
	Lef80=luminosity_R(np.asarray(Ref80),np.asarray(In2Def),\
		np.asarray(n2D_eff_2))
	LRef20=Lef20/np.asarray(Ref20)**2.0
	LRef40=Lef40/np.asarray(Ref40)**2.0
	LRef60=Lef60/np.asarray(Ref60)**2.0
	LRef80=Lef80/np.asarray(Ref80)**2.0
	print '--------------- PLOTING RESULTS -------------------'
	plt.plot(nbd,LRef20,'--o',color='red',\
		label='$(L_{	20}/R_{	20}^2)_{	ef}$')
	plt.plot(nbd,LRef40,'--o',color='blue',\
		label='$(L_{	40}/R_{	40}^2)_{	ef}$')
	plt.plot(nbd,LRef60,'--o',color='green',\
		label='$(L_{	60}/R_{	60}^2)_{	ef}$')
	plt.plot(nbd,LRef80,'--o',color='black',\
		label='$(L_{	80}/R_{	80}^2)_{	ef}$')
	# plt.xlabel('Indice de Sersic efetivo $n_{ef}$')
	plt.xlabel('$n$')
	plt.ylabel('Razao luminosidade raio: $L_i/R_i^2$')
	plt.grid()
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.title('Modelos de bojos obtidos de '+what_n[1]+'.',x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.93,x=0.70 , fontsize=13)
	fig.set_size_inches(12.0, 7.0)
	plt.savefig(path2+'Li_Ri_ef_'+what_n[0]+'.jpg')
	plt.savefig(path2+'Li_Ri_ef_'+what_n[0]+'.pdf')
	# plt.show()
	plt.clf()
	print '--------------- 40% -------------------'

# ''' --------------'''
	plt.plot(nbd,R21,'--o',color='red',label='$R_{	21n_{	BD}}$')
	plt.plot(nbd,R32,'--o',color='blue',label='$R_{	32n_{	BD}}$')
	plt.plot(nbd,R43,'--o',color='green',label='$R_{	43n_{	BD}}$')
	plt.plot(nbd,R54,'--o',color='black',label='$R_{	54n_{	BD}}$')
	plt.plot(nbd,R65,'--o',color='purple',label='$R_{	65n_{	BD}}$')
	plt.plot(nbd,R76,'--o',color='brown',label='$R_{	76n_{	BD}}$')
	plt.plot(nbd,R87,'--o',color='yellow',label='$R_{	87n_{	BD}}$')
	plt.plot(nbd,R98,'--o',color='pink',label='$R_{	98n_{	BD}}$')
	plt.plot(nbd,Ref21,'--bs',color='red',label='$R_{	21n_{	ef}}$')
	plt.plot(nbd,Ref32,'--bs',color='blue',label='$R_{	32n_{	ef}}$')
	plt.plot(nbd,Ref43,'--bs',color='green',label='$R_{	43n_{	ef}}$')
	plt.plot(nbd,Ref54,'--bs',color='black',label='$R_{	54n_{	ef}}$')
	plt.plot(nbd,Ref65,'--bs',color='purple',label='$R_{	65n_{	ef}}$')
	plt.plot(nbd,Ref76,'--bs',color='brown',label='$R_{	76n_{	ef}}$')
	plt.plot(nbd,Ref87,'--bs',color='yellow',label='$R_{	87n_{	ef}}$')
	plt.plot(nbd,Ref98,'--bs',color='pink',label='$R_{	98n_{	ef}}$')
	# plt.xlabel('Indice de Sersic efetivo $n_{	ef}$')
	plt.xlabel('$n$')
	plt.ylabel('Razao entre os raios.')
	plt.title('Modelos de bojos obtidos de'+what_n[1]+'.',x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.93,x=0.70 , fontsize=13)
	fig.set_size_inches(12.0, 7.0)
	plt.grid()
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.savefig(path2+'R_ratios_1_'+what_n[0]+'_.jpg',dpi=300)
	plt.savefig(path2+'R_ratios_1_'+what_n[0]+'_.pdf',dpi=300)
	# plt.show()
	plt.clf()
	print '--------------- 60% -------------------'
	R42=np.asarray(R40)/np.asarray(R20)
	R64=np.asarray(R60)/np.asarray(R40)
	R86=np.asarray(R80)/np.asarray(R60)
	R97=np.asarray(R90)/np.asarray(R70)
	R75=np.asarray(R70)/np.asarray(R50)
	R53=np.asarray(R50)/np.asarray(R30)
	R31=np.asarray(R30)/np.asarray(R10)
	Ref42=np.asarray(Ref40)/np.asarray(Ref20)
	Ref64=np.asarray(Ref60)/np.asarray(Ref40)
	Ref86=np.asarray(Ref80)/np.asarray(Ref60)
	Ref97=np.asarray(Ref90)/np.asarray(Ref70)
	Ref75=np.asarray(Ref70)/np.asarray(Ref50)
	Ref53=np.asarray(Ref50)/np.asarray(Ref30)
	Ref31=np.asarray(Ref30)/np.asarray(Ref10)
	plt.plot(nbd,R42,'--o',color='red',label='$R_{	42n_{	BD}}$')
	plt.plot(nbd,R64,'--o',color='blue',label='$R_{	64n_{	BD}}$')
	plt.plot(nbd,R86,'--o',color='green',label='$R_{	86n_{	BD}}$')
	plt.plot(nbd,R97,'--o',color='black',label='$R_{	97n_{	BD}}$')
	plt.plot(nbd,R75,'--o',color='purple',label='$R_{	75n_{	BD}}$')
	plt.plot(nbd,R53,'--o',color='brown',label='$R_{	53n_{	BD}}$')
	plt.plot(nbd,R31,'--o',color='yellow',label='$R_{	31n_{	BD}}$')
	plt.plot(nbd,Ref42,'--bs',color='red',label='$R_{	42n_{	ef}}$')
	plt.plot(nbd,Ref64,'--bs',color='blue',label='$R_{	64n_{	ef}}$')
	plt.plot(nbd,Ref86,'--bs',color='green',label='$R_{	86n_{	ef}}$')
	plt.plot(nbd,Ref75,'--bs',color='black',label='$R_{	75n_{	ef}}$')
	plt.plot(nbd,Ref53,'--bs',color='purple',label='$R_{	53n_{	ef}}$')
	plt.plot(nbd,Ref31,'--bs',color='brown',label='$R_{	31n_{	ef}}$')
	# plt.xlabel('Indice de Sersic efetivo $n_{	ef}$')
	plt.xlabel('$n$')
	plt.ylabel('Razao entre os raios.')
	plt.grid()
	fig.set_size_inches(12.0, 7.0)
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.title('Modelos de bojos obtidos de'+what_n[1]+'.',x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.93,x=0.70 , fontsize=13)
	plt.savefig(path2+'R_ratios_2_'+what_n[0]+'_.jpg',dpi=300)
	plt.savefig(path2+'R_ratios_2_'+what_n[0]+'_.pdf',dpi=300)
	# plt.show()
	plt.clf()
	print '--------------- 80% -------------------'
	''' Concavity relations'''
	print '------------- 	CALCULATING CONCAVITY RELATIONS -----------------'
	R=[[R21,'R21'],[R32,'R32'],[R43,'R43'],[R54,'R54'],[R65,'R65'],\
	[R76,'R76'],[R87,'R87'],[R98,'R98']]
	Ref=[Ref21,Ref32,Ref43,Ref54,Ref65,Ref76,Ref87,Ref98]
	color=['red','purple','green','black','blue','brown',\
	'yellow','gray','pink']
	# # DL=[]
	# print (nb-nbd), (nbd -n2D_eff_1) '''compare the difference between the two'''
	for j in range(1,7):
		dL=(-(R[j+1][0]-2*R[j][0]-R[j-1][0])/((R[j+1][0]-R[j][0])*(R[j][0]-R[j-1][0])))
		dLef=(abs(-(Ref[j+1]-2*Ref[j]-Ref[j-1])/((Ref[j+1]-Ref[j])*(Ref[j]-Ref[j-1]))))
		# DL.append(dL)
		plt.plot(nbd,np.log(abs(dL)),'--o',color=color[j-1],\
			label='$i='+str(j+1)+'$, $n_{	ef}$, $d^2L_{	BD}$')
		plt.plot(nbd,np.log((dLef)),'--o',color=color[j+2],\
			label='$i='+str(j+1)+'$, $n_{{	ef}}$, $d^2L_{	ef}$')
		plt.ylabel('$\ ln \left\{| d^2L_i|\\right\} $')
		# plt.xlabel('Indice de Sersic efetivo $n_{	ef}$')
		plt.xlabel('$n$')
		plt.grid()
		fig.set_size_inches(12.0, 7.0)
		plt.legend(loc='',fancybox=True, shadow=True)
		plt.title('Modelos de bojos obtidos de'+what_n[1]+'.',x=0.25)
		plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
			y=0.93,x=0.70 , fontsize=13)
		plt.savefig(path2+'nef_ln_dL_'+what_n[0]+'_i'+str(j+1)+'.pdf')
		plt.savefig(path2+'nef_ln_dL_'+what_n[0]+'_i'+str(j+1)+'.jpg')
		# plt.show()
		plt.clf()
	plt.clf()
	''' Plot all together in a same graphic.'''
	for j in range(1,7):
		dL=(-(R[j+1][0]-2*R[j][0]-R[j-1][0])/((R[j+1][0]-R[j][0])*(R[j][0]-R[j-1][0])))
		dLef=(abs(-(Ref[j+1]-2*Ref[j]-Ref[j-1])/((Ref[j+1]-Ref[j])*(Ref[j]-Ref[j-1]))))
		# DL.append(dL)
		plt.plot(nbd,np.log(abs(dL)),'--o',color=color[j-1],label='$i='+\
			str(j+1)+'$,  $n_{	ef}$, $d^2L_{	BD}$')
		plt.plot(nbd,np.log((dLef)),'--bs',color=color[j+2],label='$i='+\
			str(j+1)+'$,  $n_{{	ef}}$, $d^2L_{	ef}$')
	plt.ylabel('$\ ln \left\{| d^2L_i|\\right\} $')
	# plt.xlabel('Indices de Sersic: bojo e disco $n_{	BD}$, e \
	# efetivo $n_{	ef}$')
	# plt.xlabel('Indice de Sersic efetivo $n_{	ef}$')
	plt.xlabel('$n$')
	plt.grid()
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.title('Modelos de bojos obtidos de'+what_n[1]+'.',x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.93,x=0.70 , fontsize=13)
	fig.set_size_inches(12.0, 7.0)
	# plt.plot()
	# plt.show()
	plt.savefig(path2+'nef_ln_dL_'+what_n[0]+'_all.pdf')
	plt.savefig(path2+'nef_ln_dL_'+what_n[0]+'_all.jpg',dpi=300)
	plt.clf()
	print '--------------- 100% -------------------'
	# pgBar.stop()
	print '--------------- DONE -------------------'
	return 

def compare_mfmtk_parameters(path,path2,var1,var2,var3,plot_config):
	'''
		This function recieves the same morphometric paramaters, 
		calculated by morfometryka, in two different galaxies: 
		one with a bulge+disk component, and the other, a bulge, created 
		with the effective Sersic of the bulge+disk. 

		It saves many graphic relations.
	'''
	# var1=['n1D',n1D,'$n_{1D}$']
	# var2=['C1',C1,'$C_1$']
	# var3=['C2',C2,'$C_2$']
	fig = plt.gcf()
	plt.plot(var1[1],var2[1],'--o',color='red',label='$'+var2[2]+'$')
	plt.plot(var1[1],var3[1],'--o',color='blue',label='$'+var3[2]+'$')
	plt.grid()
	plt.legend(loc='',fancybox=True, shadow=True)
	plt.xlabel('$'+var1[2]+'$')
	plt.ylabel('$'+var2[2]+'$  e  $'+var3[2]+'$')
	plt.title('$'\
			  +str(plot_config[0][2])+'='+str(plot_config[0][1][0])+'$, $'\
			  +str(plot_config[1][2])+'='+str(plot_config[1][1][0])+'$, $'\
			  +str(plot_config[2][2])+'='+str(plot_config[2][1][0])+'$, $'\
			  +str(plot_config[3][2])+'='+str(plot_config[3][1][0])+'$, $'\
			  +str(plot_config[4][2])+'='+str(plot_config[4][1][0])+'$',\
			  x=0.25)
	plt.suptitle(r'Medidas Morfometryka G2 (Ferrari, F., 2015)', \
		y=0.92,x=0.70 , fontsize=13)
	fig.set_size_inches(12.0, 7.0)
	plt.savefig(path2+'comparisons/'+var1[0]+'_'+var2[0]+var3[0]+\
		'.jpg',dpi=300)
	plt.savefig(path2+'comparisons/'+var1[0]+'_'+var2[0]+var3[0]+\
		'.pdf',dpi=300)
	print 'Plot '+var2[0]+' e '+var3[0]+' versus '+var1[0]+' done!'
	# plt.show()
	plt.clf()


	

##############################################################################
## MORFOMETRYKA UTILS - Ferrari, F. et al, 2015
##############################################################################
def polarim(image, origin=None, log=False ):
      """Reprojects a 2D numpy array ("image") into a polar coordinate system.
      "origin" is a tuple of (x0, y0) and defaults to the center of the image.
      http://stackoverflow.com/questions/3798333/image-information-along-a-polar-coordinate-system
      refactored by FF, 2013-2014 (see transpolar.py)
      """

      if origin == None:
         origin = np.array(image.shape)/2.


      def cart2polar(x, y):
         r = np.sqrt(x**2 + y**2)
         theta = np.arctan2(y, x)
         return r, theta

      def polar2cart(r, theta):
         x = r * np.cos(theta)
         y = r * np.sin(theta)
         return x, y

      def cart2logpolar(x, y, M=1):
         alpha  = 0.01
         r = np.sqrt(x**2 + y**2)
         rho = M * np.log(r + alpha)
         theta = np.arctan2(y, x)
         return rho, theta

      def logpolar2cart(rho, theta, M=1):
         x = np.exp(rho/M) * np.cos(theta)
         y = np.exp(rho/M) * np.sin(theta)
         return x, y

        
      ny, nx = image.shape 
      if origin is None:
         x0, y0 = (nx // 2, ny // 2)
         origin = (x0,y0)
      else:
         x0, y0 = origin

      # Determine that the min and max r and theta coords will be...
      x, y = np.meshgrid( np.arange(nx) - x0, np.arange(ny) - y0 )

      r, theta = cart2polar(x, y)

      # Make a regular (in polar space) grid based on the min and max r & theta
      r_i     = np.linspace(r.min(),     r.max(),     nx )
      theta_i = np.linspace(theta.min(), theta.max(), ny)
      theta_grid, r_grid = np.meshgrid(theta_i, r_i)

      # Project the r and theta grid back into pixel coordinates
      xi, yi = polar2cart(r_grid, theta_grid)
      xi += origin[0] # We need to shift the origin back to 
      yi += origin[1] # back to the lower-left corner...
      xi, yi = xi.flatten(), yi.flatten()
      coords = np.vstack((xi, yi)) # (map_coordinates requires a 2xn array)
      
      zi = nd.map_coordinates(image, coords, order=1)
      galpolar      = zi.reshape((nx,ny))

      # self.r_polar  = r_i
      # self.theta_polar = theta_i


      # return galpolar, r, theta, theta_grid,r_grid
      return galpolar

def grad(polar_data):
	dx,dy = savitzky_golay_2d(polar_data, 11,3,'both')
	ort = np.arctan2(dy, dx)
	ortn =  (ort+np.pi)%(np.pi)
	# plot(std(ortn, axis=1))
	return np.std(ortn, axis=1)


def lum_T(In,Rn,n):
	return ((2*np.pi*In*Rn**2.0*n*np.exp(bn(n)))/(bn(n)**(2.0*n)))*gamma(2*n)
def Hen(In,Rn,n):
	# return -lum_T(In,Rn,n)*(np.log(In)+bn(n)+gamma(2*n+1)/(bn(n)*gamma(2*n))	)
	return -(np.log(In)+bn(n)+gamma(2*n+1)/(bn(n)*gamma(2*n))	)