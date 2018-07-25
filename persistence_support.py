# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0


import os,sys,glob,time,collections,gc
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import cPickle as pickle
import dimarray as da
from scipy.optimize import curve_fit
from lmfit import  Model
import pandas as pd

# ---------- Counter Conversion ------------------------
def counter_to_list(counter):
	tmp=[]
	lengths=counter.keys()
	if 0 in lengths:
		lengths.remove(0)
	if len(lengths)>2:
		for key in lengths:
			for i in range(counter[key]):
				tmp.append(key)
		tmp=np.array(tmp)
		return -tmp[tmp<0],tmp[tmp>0]
	else:
		return [],[]


def counter_to_pers(counter,state):
	# to distribution
	pers_tmp=np.array(range(1,max([state*key+1 for key in counter.keys()])))
	count,pers=[],[]
	for pp,i in zip(pers_tmp,range(len(pers_tmp))):
		if pp*state in counter.keys():
			if abs(pp*state)!=0:
				count.append(counter[pp*state])
				pers.append(pp)
	count=np.array(count)
	pers=np.array(pers)
	return count,pers


# ---------- Fitting -----------------------------------
def double_exp(x,a1,b1,b2,thr):
	x=np.array(x)
	y=x.copy()*np.nan
	y[x<=thr]	=	a1*np.exp(-b1*x[x<=thr])
	y[x>thr]	=	a1*np.exp((b2-b1)*thr)*np.exp(-b2*(x[x>thr]))
	return y

def single_exp(x,a1,b1):
	x=np.array(x)
	return a1*np.exp(-b1*x)

def two_exp(x,a1,b1,a2,b2):
	x=np.array(x)
	y=x.copy()*np.nan
	y	=	a1*np.exp(-b1*x)+a2*np.exp(-b2*x)
	return y

def all_fits(count,pers,plot=False,subax=None,do_two_expo=False):
	tmp={}

	# single
	try:
		popt, pcov = curve_fit(single_exp,pers[2::],count[2::])
		model=Model(single_exp)
		model.set_param_hint('a1', value=popt[0],vary=False)
		model.set_param_hint('b1', value=popt[1],vary=False)
		result = model.fit(count[2::], x=pers[2::])
		tmp['single_exp']={'bic':result.bic,'params':result.best_values,'best_fit':result.best_fit}
	except Exception,e:
		tmp['single_exp']={'bic':None,'params':{'a1':np.nan,'b1':np.nan},'best_fit':None}


	# double
	try:
		popt_, pcov = curve_fit(double_exp,pers[2::],count[2::],p0=[popt[0],popt[1],popt[1],7.],bounds=([0,0,0,5.],[np.inf,np.inf,np.inf,14.]))
		doubleM=Model(double_exp)
		doubleM.set_param_hint('a1', value=popt_[0],vary=False)
		doubleM.set_param_hint('b1', value=popt_[1],vary=False)
		doubleM.set_param_hint('b2', value=popt_[2],vary=False)
		doubleM.set_param_hint('thr', value=popt_[3],vary=False)
		result = doubleM.fit(count[2::], x=pers[2::])
		tmp['double_exp']={'bic':result.bic,'params':result.best_values,'best_fit':result.best_fit}
	except Exception,e:
		tmp['double_exp']={'bic':None,'params':{'a1':np.nan,'b1':np.nan,'b2':np.nan,'thr':np.nan},'best_fit':None}

	if do_two_expo:
		try:
			popt_, pcov = curve_fit(two_exp,pers[2::],count[2::],p0=[popt[0]*0.1,popt[1],popt[0]*10,popt[1]],bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf]))
			doubleM=Model(two_exp)
			doubleM.set_param_hint('a1', value=popt_[0],vary=False)
			doubleM.set_param_hint('b1', value=popt_[1],vary=False)
			doubleM.set_param_hint('a2', value=popt_[2],vary=False)
			doubleM.set_param_hint('b2', value=popt_[3],vary=False)
			result = doubleM.fit(count[2::], x=pers[2::])
			tmp['two_exp']={'bic':result.bic,'params':result.best_values,'best_fit':result.best_fit}
		except Exception,e:
			tmp['two_exp']={'bic':None,'params':{'a1':np.nan,'b1':np.nan,'a2':np.nan,'b2':np.nan},'best_fit':None}

	if plot:
		subax.plot(pers[2::],count[2::])
		subax.plot(pers[2::],tmp['single_exp']['best_fit'],label='single '+str(round(tmp['single_exp']['bic'],2)))
		subax.plot(pers[2::],tmp['double_exp']['best_fit'],label='double '+str(round(tmp['double_exp']['bic'],2)))
		subax.plot(pers[2::],tmp['two_exp']['best_fit'],label='two '+str(round(tmp['two_exp']['bic'],2)))
		subax.set_yscale('log')
		subax.set_xlim((0,40))
		subax.set_ylim((100,count[2]))
		subax.set_title(region)

	return tmp

# fig = plt.figure(figsize=(9,6))
# ax_big=fig.add_axes([0,0,1,1])
# all_fits(region_dict['NAS']['DJF']['cold']['count'],region_dict['NAS']['DJF']['cold']['pers'],plot=True,subax=ax_big)
# plt.savefig('plots/NAS.png')
