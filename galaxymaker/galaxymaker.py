# !/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
## This is a main script that runs routines to create synthetic galaxies    ##
## 			Created by Geferson Lucatelli 			    ##
##			lucatelli_asph_r2@virgophysics 	 		    ##
##					  			            ##
##		(original version 03.16.2016) (lastest version 05.23.2016)  ##
##############################################################################

##############################################################################
##IMPORT SCRIPTS
##############################################################################
from galaxymakerlib import*
import os
import requests
import sys
import time
from multiprocessing import Process
from mfmtkutils import *
"""Create synthetic galaxies, including in ellipticals and spirals\
components like bulges, disks and bars.\
"""
"""The function create_galaxies() calls a specified profile, to create a\
   synthetic galaxie. \
   
   Returns an image array\

   Some components that can be used to control the intensity of the\
   profile in the galaxy are:\
		
		ec	exponential contribution to the final galaxy intensity output\
		
		buc	bulge contribution to the final galaxy intensity output\
		
		dc 	disk contribution to the final galaxy intensity output\
		
		barc	bar contribution to the final galaxy intensity output\
		
		nsc	general sersic profile contribution to the final galaxy 
		intensity output
		
		sct	tangent spiral contribution to the final galaxy intensity \
		output
		
		scth	hyperbolic tangent spiral contribution to the final\
		galaxy intensity output

		gal_center_ass 		include a profile in any other center. \
		It is useful to create assimetric galaxies.



		It can be created a simple synthetic galaxie in a more easy \
		way simply calling a function.\
"""

def create_galaxies(bd,nn,rn,Re,L_T,c,qb,Qd,SPIRAL,path):
	plt.gray()
	psf=pf.getdata('psf_final.fits')
	nsc=1.0
	ec=1.0
	spo=2.0
	PA=0.0
	number_name=1
	gal_center=(0.0,0.0) #for all components.
	''' Make the galaxies.'''
	if SPIRAL is True:
		''' Spiral Quantities'''
		# N=[0.77,1.0,1.24,1.53,1.78,2.05,2.32,2.57,2.85,3.10]
		PHI0=[0.69,0.98]
		# N=[1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8]
		N=[1.65,1.84]
		# PHI0=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
		AA=1.3
		p=1.0
		k=2.0
		LSP=[]  #to store the values of the total spiral luminosity
		for phi0 in PHI0:
			for NN in N:
				for re in Re:
					for BD in bd:
						start_time = time.time()
						# BT=BD*((1+BD)**(-1.0))
						# Ie,L_BT,L_DT,L_T,BDd=bulge_disk_component(n,rn,In,re,BT)
						L_BT,L_DT,BT,In,Ie=bulge_disk(BD,n,rn,re,L_T)
						elliptical,sersic,disk=ell_gal(n,In,rn,Ie,re,ec,nsc,\
							qb,qd,c,PA,gal_center)
						spiral=tan_spiral_profile(k,spo,AA,NN,phi0,q,PA,\
							gal_center,re,Ie)
						spiral=convolve(spiral,psf)
						sersic=convolve(sersic,psf)
						disk=convolve(disk,psf)
						noise=1.0*np.random.rand(255,255)
						# lT=spiral.sum()+disk.sum()+sersic.sum()
						# ls=spiral.sum()/lT
						# ld=disk.sum()/lT
						# lb=sersic.sum()/lT
						# LSP.append(ls)
						gal=sersic+disk+spiral+noise
						plot_save_gal(gal,sersic,disk,spiral,n,In,rn,Ie,re,\
							q,c,PA,BT,BD,L_BT/L_T,L_DT/L_T,L_T/L_T,nsc,ec,\
							k,spo,AA,NN,phi0,number_name,path,SPIRAL)
						plot_grad(gal,sersic,disk,spiral,number_name,path,\
							SPIRAL=True)
						save_fits(gal,number_name,None,path)
						save_fits(polarim(gal)[:150],number_name-1,\
							'gal_polar_'+str(number_name),path)
						w_values('gal',(n,In,rn,Ie,re,q,c,PA,BT,BD,L_BT/L_T,\
							L_DT/L_T,L_T/L_T,nsc,ec,C1(n),\
						C2(n),Crn(n),RP(n,rn),k,spo,AA,NN,phi0),path)
						print 'Gal ', number_name, ' created', 'Time= ', \
						"--- %s seconds ---" % (time.time() - start_time)
						number_name=number_name+1
		# return LSP
	else:
		SPIRAL=False
		"""Set to unity the spiral input parameters, it is for convenience."""
		NN=1.0
		phi0=1.0
		AA=1.0
		p=1.0
		k=1.0
		for qd in Qd:
			for n in nn[0]:
				for re in Re:
					for BD in bd:
						start_time = time.time()
						
						L_BT,L_DT,BT,In,Ie,mu_n,mean_mu_n,mu_e,mean_mu_e=\
						bulge_disk(BD,n,rn,re,L_T)
						
						elliptical,sersic,disk=ell_gal(n,In,rn,Ie,re,ec,\
							nsc,qb,qd,c,PA,gal_center)
						
						spiral=0.0
						sersic=convolve(sersic,psf)
						disk=convolve(disk,psf)
						noise=1.0*np.random.rand(255,255)
						gal=sersic+disk+noise
						
						plot_save_gal(gal,sersic,disk,spiral,n,nn,In,rn,Ie,\
							re,qb,qd,c,PA,BT,BD,L_BT/L_T,L_DT/L_T,L_T/L_T,\
							nsc,ec,k,p,AA,NN,phi0,number_name,path,\
							SPIRAL,DISK=True)
						save_fits(gal,number_name,None,path)
						
						save_fits(polarim(gal)[:150],number_name,\
							'gal_polar_'+str(number_name),path)
						
						plot_grad(gal,sersic,disk,spiral,number_name,\
							path,SPIRAL)
						
						w_values('gal',('gal_'+str(number_name),n,In,rn,\
							Ie,re,mu_n,mean_mu_n,mu_e,mean_mu_e,qb,qd,c,\
							PA,BT,BD,L_BT/L_T,L_DT/L_T,L_T/L_T,nsc,ec,\
							C1f(n),C2f(n),Crn(n),RPe(n,rn),k,p,AA,NN,phi0)\
							,path)
						
						print 'Gal ', number_name, ' created', 'Time= ', \
						"--- %s seconds ---" % (time.time() - start_time)
						number_name=number_name+1
	"""Return the total number of galaxies created."""
	return len(nn[0])*len(bd)

def create_gal_eff_n(nfit,path):
	"""Create a set of galaxies with the effective Sersic index calculated\
	with MORFOMETRYKA. 
	"""
	DISK=False
	plt.gray()
	psf=pf.getdata('psf_final.fits')
	''' Sersic Quantitites'''
	rn=10.0
	L_T=3.0e5
	re=0.0
	nsc=1.0
	ec=1.0
	spo=2.0 
	Ie=0.0 #No disk.
	BD=1.0 #Set to 1 to for convenience
	BT=1.0
	L_BT=L_T
	L_DT=0.0
	''' General Ellipsis Quantities'''
	qb=1.0
	qd=1.0
	c=0.0
	PA=0.0
	number_name=1
	gal_center=(0.0,0.0) #for all components.
	SPIRAL=False
	NN=1.0
	phi0=1.0
	AA=1.0
	p=1.0
	k=1.0
	for n in nfit[0]:
		start_time = time.time()
		In=Inn(n[0],rn,L_BT)
		
		elliptical,sersic,disk=ell_gal(n[0],In,rn,Ie,re,ec,nsc,\
			qb,qd,c,PA,gal_center)
		
		spiral=0.0
		disk[:]=0.0
		mu_n=mean_mu_n=mu_e=mean_mu_e=0
		# print (sersic+disk).max()
		
		sersic=convolve(sersic,psf)
		noise=1.0*np.random.rand(255,255)
		gal=sersic+noise
		
		plot_save_gal(gal,sersic,disk,spiral,n[0],[nfit[0],nfit[1]],\
			In,rn,Ie,re,qb,qd,c,PA,BT,BD,L_BT/L_T,\
			L_DT/L_T,L_T/L_T,nsc,ec,k,p,AA,NN,phi0,number_name,path,\
			SPIRAL=False,DISK=False)
		
		save_fits(gal,number_name,None,path)
		save_fits(polarim(gal)[:150],number_name,'gal_polar_'+\
			str(number_name),path)
		
		plot_grad(gal,sersic,disk,spiral,number_name,path,SPIRAL=False)

		w_values('gal',('gal_'+str(number_name),n[0],In,rn,Ie,re,mu_n,\
			mean_mu_n,mu_e,mean_mu_e,qb,qd,c,PA,BT,BD,L_BT/L_T,L_DT/L_T,\
			L_T/L_T,nsc,ec,C1f(n[0]),C2f(n[0]),Crn(n[0]),RPe(n[0],rn),k,\
			p,AA,NN,phi0),path)
		
		print 'Gal ', number_name, ' created', 'Time= ',\
		"--- %s seconds ---" % (time.time() - start_time)
		number_name=number_name+1
	"""Return the total number of galaxies."""
	return len(nfit[0])

"""Funcions to call the morfometrykaG2 code.
"""
def mfmtk1(path,nber_galaxies):
	for i in range(1,10):
		start_time = time.time()
		os.system("python /home/lucatelli/galaxy_maker/morfometrykaG2.py "+path+"fits/gal_0"+\
			str(i)+".fits psf_final.fits noshow")
		print 'Gal ', str(i),' done.'
		print"--- %s seconds ---" % (time.time() - start_time)
		pass
def mfmtk2(path,nber_galaxies):
	for i in range(10,nber_galaxies+1):
		start_time = time.time()
		os.system("python /home/lucatelli/galaxy_maker/morfometrykaG2.py "+path+"fits/gal_"+\
			str(i)+".fits psf_final.fits noshow")
		print 'Gal ', str(i),' done.'
		print"--- %s seconds ---" % (time.time() - start_time)
		pass
def runInParallel(*fns):
	proc = []
	for fn in fns:
		p = Process(target=fn)
		p.start()
		proc.append(p)
	for p in proc:
		p.join()

def callPlots(mfmtk_cat,glm_cat,path,varss,cruzed_values=True):
	'''Plot many quantities with some configuration.
	   
	   varss example:
	   
	   		varrs[0] correspond to the fundamental variable, and
	   		varrs[1] correspond to the constante one.

	   		varss=[	[['BD'],['\\xi_{BD}']]	,['n',['n'],'n']]
	   		varss=[	[['n'],['n']]	,['BD',['BD'],'\\xi_{BD}']]
	'''
	gml_par=varss[0][0]
	mfmtk_pars=[['C1'],['C2'],['H'],['G'],['A1'],['A3'],['S1'],['S3'],\
				['M20'],['RnFit1D'],['RnFit2D'],['InFit1D'],['InFit2D'],\
				['R10'],['R20'],['R30'],['R40'],['R50'],['R60'],['R70'],\
				['R80'],['R90'],['psi'],['sigma_psi'],['nFit1D'],\
				['nFit2D']]
	labels=[varss[0][1][0],\
	        ['C_1','C_2','H','G','A_1','A_3','S_1','S_3','M_{20}',\
	        'R_{n_{1D}}','R_{n_{2D}}','I_{n_{1D}}','I_{n_{2D}}',\
	        'R_{10}','R_{20}','R_{30}','R_{40}','R_{50}','R_{60}',\
	        'R_{70}','R_{80}','R_{90}','\psi','\sigma_\psi',\
	        'n_{1D}','n_{2D}']]
	plot_config=[
				['qb',glm_cat.param_selection(['qb'],\
					column_dict_glm)[0],'q_b'],\
				['qd',glm_cat.param_selection(['qd'],\
					column_dict_glm)[0],'q_d'],\
				[varss[1][0] ,glm_cat.param_selection(varss[1][1],\
					column_dict_glm)[0],varss[1][2]],\
				['Rn',glm_cat.param_selection(['Rn'],\
					column_dict_glm)[0],'R_n'],\
				['RdRn',glm_cat.param_selection(['Rd'],\
					column_dict_glm)[0]/glm_cat.param_selection(['Rn'],\
				 column_dict_glm)[0],'R_d/R_n']]
	for i in range(len(mfmtk_pars)):
		plots(glm_cat.param_selection(gml_par,column_dict_glm),\
			  mfmtk_cat.param_selection(mfmtk_pars[i],column_dict),\
			  [labels[0],labels[1][i]],plot_config,path)
	'''If we want to study how the measuerd parameters correlates each other, 
	lets plot they values in pairs.
	'''
	if cruzed_values is True:
		print '-------------- PLOTTING CROSSED PARAMETERS --------------'
		s=1
		for j in range(len(mfmtk_pars)):
			for k in range(s,len(mfmtk_pars)):#do not repeat the plots!
				plots(mfmtk_cat.param_selection(mfmtk_pars[j],column_dict),\
					mfmtk_cat.param_selection(mfmtk_pars[k],column_dict),\
					[labels[1][j],labels[1][k]],plot_config,path)
				# print [labels[1][j],labels[1][k]]
			s=j+2
	return 'Done!'

start_time0 = time.time()
path='/home/lucatelli/galaxy_maker/data_galaxies/convolved_noise_2/n_var/BD_8/'
path2=path+'n2eff/'
path3=path+'n1eff/'
nn=[np.arange(1.80,7.00,2.0),'n']
bd=[3.0]
rn=10
L_T=3.0e5
Re=[3.0*rn]
c=0
qb=1.0
Qd=[1.0]
SPIRAL=False
#########################################################################
######################    		STEP #1			  #######################
#########################################################################
print '-------------- PREPARING TO CREATE THE GALAXIES --------------'
'''Create a set of synthetic galaxies.'''
nber_galaxies=create_galaxies(bd,nn,rn,Re,L_T,c,qb,Qd,SPIRAL,path)
print '-------------- PREPARING TO CALL THE MORFOMETRYKA CODE --------------'
'''Call the morfometryka code.'''
mfmtk1(path,nber_galaxies)
mfmtk2(path,nber_galaxies)

print '-------- READING GALAXYMAKER AND MORFOMETRYKA OUTPUT FILES --------'

'''Load the stored values of parameters. 
   They are the galaxymaker inputs and the morfometryka outputs.
'''
os.chdir(path+'fits/')
os.system("cat *.mfmtk | sed '/#/d' | sed 's/\\\//g'  > gal.mfmtk")
mfmtk_cat= catalog(path=path+'fits/gal.mfmtk')
glm_cat= catalog(path=path+'gal_values.dat')
os.chdir('/home/lucatelli/galaxy_maker/')
print '-------------- PREPARING TO PLOT --------------'
varss=[	[['n'],['n']]	,['BD',['BD'],'\\xi_{BD}']]
callPlots(mfmtk_cat,glm_cat,path,varss)
#########################################################################
######################    		STEP #2			  #######################
#########################################################################
'''
	Creating The New Set Of Galaxies With The n2D \
	Sersic Index Calculated In Step #1
'''
print '-------------- LOADING EFFECTIVE SERSIC INDEX --------------'
nfit2D=[mfmtk_cat.param_selection(['nFit2D'],column_dict),'n2D']
what_n=['n','$n$']
print '-- PREPARING TO MAKE NEW GALAXIES WITH THE EFFECTIVE SERSIC INDEX --'
nber_galaxies=create_gal_eff_n(nfit2D,path2)
print '-------------- PREPARING TO CALL THE MORFOMETRYKA CODE --------------'
mfmtk1(path2,nber_galaxies)
mfmtk2(path2,nber_galaxies)
print '--------- READING GALAXYMAKER AND MORFOMETRYKA OUTPUT FILES ---------'
'''Load the stored values of parameters. 
   They are the galaxymaker inputs and the morfometryka outputs.
'''
os.chdir(path2+'fits/')
os.system("cat *.mfmtk | sed '/#/d' | sed 's/\\\//g'  > gal.mfmtk")
os.chdir('/home/lucatelli/galaxy_maker/')
mfmtk_cat2= catalog(path=path2+'fits/gal.mfmtk')
glm_cat2= catalog(path=path2+'gal_values.dat')
varss=[	[['n'],['n']]	,['BT',['BT'],'\\xi_{BT}']]
print '-------------- PREPARING TO PLOT --------------'
callPlots(mfmtk_cat2,glm_cat2,path2,varss)
compare_n_R(path,path2,what_n,mfmtk_cat,glm_cat,mfmtk_cat2,glm_cat2)

"""Compare morphometric parameters."""
mfmtk_pars=[['C1'],['C2'],['H'],['G'],['A1'],['A3'],['S1'],['S3'],\
			['M20'],['RnFit1D'],['RnFit2D'],['InFit1D'],['InFit2D'],\
			['R10'],['R20'],['R30'],['R40'],['R50'],['R60'],['R70'],\
			['R80'],['R90'],['psi'],['sigma_psi'],['nFit1D'],['nFit2D']]
var1=['n',glm_cat.param_selection(['n'][0],column_dict_glm),'n']
labels=[['n'],\
		['C_1','C_2','H','G','A_1','A_3','S_1','S_3','M_{20}','R_{n_{1D}}',\
		'R_{n_{2D}}','I_{n_{1D}}','I_{n_{2D}}','R_{10}','R_{20}','R_{30}',\
		'R_{40}','R_{50}','R_{60}','R_{70}','R_{80}','R_{90}','\psi',\
		'\sigma_\psi','n_{1D}','n_{2D}']]
labels_ef=[['n'],\
		['C_{1_{ef}}','C_{2_{ef}}','H_{ef}','G_{ef}','A_{1_{ef}}',\
		'A_{3_{ef}}','S_{1_{ef}}','S_{3_{ef}}','M_{20_{ef}}',\
		'R_{n_{1Def}}','R_{n_{2Def}}','I_{n_{1Def}}','I_{n_{2Def}}',\
		'R_{10_{ef}}','R_{20_{ef}}','R_{30_{ef}}','R_{40_{ef}}',\
		'R_{50_{ef}}','R_{60_{ef}}','R_{70_{ef}}','R_{80_{ef}}',\
		'R_{90_{ef}}','\psi_{ef}','\sigma_{\psi_{ef}}','n_{1D_{ef}}',\
		'n_{2D_{ef}}']]
plot_config=[
			 ['qb',glm_cat.param_selection(['qb'],\
			 	column_dict_glm)[0],'q_b'],\
			 ['qd',glm_cat.param_selection(['qd'],\
			 	column_dict_glm)[0],'q_d'],\
			 ['BD',glm_cat.param_selection(['BD'],\
			 	column_dict_glm)[0],'\\xi_{BD}'],\
			 ['Rn',glm_cat.param_selection(['Rn'],\
			 	column_dict_glm)[0],'R_n'],\
			 ['RdRn',glm_cat.param_selection(['Rd'],\
			 	column_dict_glm)[0]/glm_cat.param_selection(['Rn'],\
			 	column_dict_glm)[0],'R_d/R_n']]
for i in range(len(mfmtk_pars)):
	var2=[mfmtk_pars[i][0],mfmtk_cat.param_selection(mfmtk_pars[i],\
		column_dict),labels[1][i]]
	var3=[mfmtk_pars[i][0],mfmtk_cat2.param_selection(mfmtk_pars[i],\
		column_dict),labels_ef[1][i]]
	compare_mfmtk_parameters(path,path2,var1,var2,var3,plot_config)


'''
	Creating The New Set Of Galaxies With The n1D Sersic Index Calculated In Step #1
'''
print '------------ LOADING EFFECTIVE SERSIC INDEX ------------'
nfit1D=[mfmtk_cat.param_selection(['nFit1D'],column_dict),'n1D']
print '-- PREPARING TO MAKE NEW GALAXIES WITH THE EFFECTIVE SERSIC INDEX --'
nber_galaxies=create_gal_eff_n(nfit1D,path3)
print '-------------- PREPARING TO CALL THE MORFOMETRYKA CODE --------------'
mfmtk1(path3,nber_galaxies)
mfmtk2(path3,nber_galaxies)
print '--------- READING GALAXYMAKER AND MORFOMETRYKA OUTPUT FILES ---------'
what_n=['n1D','$n_{	1D}$']
os.chdir(path3+'fits/')
os.system("cat *.mfmtk | sed '/#/d' | sed 's/\\\//g'  > gal.mfmtk")
os.chdir('/home/lucatelli/galaxy_maker/')
mfmtk_cat3= catalog(path=path3+'fits/gal.mfmtk')
glm_cat3= catalog(path=path3+'gal_values.dat')
varss=[	[['n'],['n']]	,['BT',['BT'],'\\xi_{BT}']]
print '-------------- PREPARING TO PLOT --------------'
callPlots(mfmtk_cat3,glm_cat3,path3,varss)
compare_n_R(path,path3,what_n,mfmtk_cat,glm_cat,mfmtk_cat3,glm_cat3)

"""Compare morphometric parameters."""
mfmtk_pars=[['C1'],['C2'],['H'],['G'],['A1'],['A3'],['S1'],['S3'],\
			['M20'],['RnFit1D'],['RnFit2D'],['InFit1D'],['InFit2D'],\
			['R10'],['R20'],['R30'],['R40'],['R50'],['R60'],['R70'],\
			['R80'],['R90'],['psi'],['sigma_psi'],['nFit1D'],['nFit2D']]
var1=['n',glm_cat.param_selection(['n'][0],column_dict_glm),'n']
labels=[['n'],\
		['C_1','C_2','H','G','A_1','A_3','S_1','S_3','M_{20}',\
		'R_{n_{1D}}','R_{n_{2D}}','I_{n_{1D}}','I_{n_{2D}}',\
		'R_{10}','R_{20}','R_{30}','R_{40}','R_{50}','R_{60}',\
		'R_{70}','R_{80}','R_{90}','\psi',\
		'\sigma_\psi','n_{1D}','n_{2D}']]
labels_ef=[['n'],\
		['C_{1_{ef}}','C_{2_{ef}}','H_{ef}','G_{ef}','A_{1_{ef}}',\
		'A_{3_{ef}}','S_{1_{ef}}','S_{3_{ef}}','M_{20_{ef}}',\
		'R_{n_{1Def}}','R_{n_{2Def}}','I_{n_{1Def}}','I_{n_{2Def}}',\
		'R_{10_{ef}}','R_{20_{ef}}','R_{30_{ef}}','R_{40_{ef}}',\
		'R_{50_{ef}}','R_{60_{ef}}','R_{70_{ef}}','R_{80_{ef}}',\
		'R_{90_{ef}}','\psi_{ef}','\sigma_{\psi_{ef}}','n_{1D_{ef}}',\
		'n_{2D_{ef}}']]
plot_config=[
			 ['qb',glm_cat.param_selection(['qb'],\
			 	column_dict_glm)[0],'q_b'],\
			 ['qd',glm_cat.param_selection(['qd'],\
			 	column_dict_glm)[0],'q_d'],\
			 ['BD',glm_cat.param_selection(['BD'],\
			 	column_dict_glm)[0],'\\xi_{BD}'],\
			 ['Rn',glm_cat.param_selection(['Rn'],\
			 	column_dict_glm)[0],'R_n'],\
			 ['RdRn',glm_cat.param_selection(['Rd'],\
			 	column_dict_glm)[0]/glm_cat.param_selection(['Rn'],\
			 	column_dict_glm)[0],'R_d/R_n']]
for i in range(len(mfmtk_pars)):
	var2=[mfmtk_pars[i][0],mfmtk_cat.param_selection(mfmtk_pars[i],\
		column_dict),labels[1][i]]
	var3=[mfmtk_pars[i][0],mfmtk_cat3.param_selection(mfmtk_pars[i],\
		column_dict),labels_ef[1][i]]
	compare_mfmtk_parameters(path,path3,var1,var2,var3,plot_config)

print '----------------------------------------------------------------------'
print '---------------------- ALL CALCULATIONS FINISHED ---------------------'
print"--- %s minutes ---" % ((time.time() - start_time0)/60.0)
print '-------- The galaxy_maker code,\
	  by Geferson Lucatelli, Fabricio Ferrari and Leonardo Ferreira ---------'
print '----------------------------------------------------------------------'