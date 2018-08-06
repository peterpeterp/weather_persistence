
# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

import os,sys,glob,time,collections
import numpy as np
from netCDF4 import Dataset,num2date
import random as random


def period_identifier(ind):
	"""
	Straight foreward period identifier. Missing values cut periods.

	Parameters
	----------
		ind: np.array
			array containing -1 and 1 corresponding to two state

	Returns
	--------
		pers: np.array
			array of the same length as `ind` containing the length of identified periods. Periods of state -1 are have negative lengths. The period length of a period is placed in the center of the period. All other values are 0.
	"""

	pers=ind.copy()*0.0

	state,count=ind[0],1
	for i in range(1,len(ind)):
		if ind[i]!=state:
			pers[i-count/2-1]=state*count
			count=0
			state*=-1
			if np.isfinite(ind[i])==False:
				pers[i]=np.nan
				count=0
		if ind[i]==state:
			count+=1

	# still an issue with last period??
	if state==1:	pers[i-count/2]=state*count
	if state==-1:	pers[i-count/2]=state*count

	return(pers)


def optimized_period_identifier(ind):
	"""
	This function identifies persistent periods using collections. It isn't a straight foreward implementation but runs faster than :meth:`period_identifier`

	Parameters
	----------
		ind: np.array
			array containing -1 and 1 corresponding to two state

	Returns
	--------
		pers: np.array
			array of the same length as `ind` containing the length of identified periods. Periods of state -1 are have negative lengths. The period length of a period is placed in the center of the period. All other values are 0.
	"""
	pers=ind.copy()*0

	ind[ind<0]=0

	cuts=list(np.where(np.isfinite(ind)==False)[0])
	cuts.append(len(ind))
	cut_start=0

	if len(np.where(ind==0)[0])==len(ind):
		return pers

	for cut_stop in cuts:
		if cut_start==cut_stop:
			cut_start=cut_stop+1
		else:
			ind_cut=ind[cut_start:cut_stop]
			pers_cut=ind_cut.copy()*0

			su=np.cumsum(ind_cut)
			counter=collections.Counter(su)

			index=0
			for count,val in zip(counter.values(),counter.keys()):
				index+=count
				if count>1:
					pers_cut[index-(count-1)/2-1]=-1*(count-1)
			# correct start
			if len(ind_cut)==1:	pers_cut[0]=-1
			else:
				if ind_cut[0]==0 and ind_cut[1]==1:	pers_cut[0]=-1
				if ind_cut[0]==0 and ind_cut[1]==0:	pers_cut[np.where(pers_cut<0)[0][0]]-=1

			ind_cut=-ind_cut+1
			su=ind_cut.copy()*0 + np.nan_to_num(ind_cut).cumsum()
			counter=collections.Counter(su)

			index=0
			for count,val in zip(counter.values(),counter.keys()):
				index+=count
				if count>1:
					pers_cut[index-(count-1)/2-1]=count-1
			# correct start
			if len(ind_cut)==1:	pers_cut[0]=1
			else:
				if ind_cut[0]==0 and ind_cut[1]==1:	pers_cut[0]=1
				if ind_cut[0]==0 and ind_cut[1]==0:	pers_cut[np.where(pers_cut>0)[0][0]]+=1

			pers[cut_start:cut_stop]=pers_cut
			cut_start=cut_stop+1

	return(pers)

def test_persistence(N):
	ind=np.random.random(N)
	ind[-2]=np.nan
	ind[ind<0.5]=-1
	ind[ind>=0.5]=1
	ind=np.array(ind,'f')
	print(ind[0:100])

	start_time = time.time()
	print(period_identifier(ind)[0:100])
	print("--- basic_and_understandable %s seconds ---" % (time.time() - start_time))

	start_time = time.time()
	print(optimized_period_identifier(ind)[0:100])
	print("--- optimized_period_identifier %s seconds ---" % (time.time() - start_time))

#test_persistence(100)

def get_persistence(state_file,out_file,seasons={'MAM':{'months':[3,4,5],'index':0}, 'JJA':{'months':[6,7,8],'index':1}, 'SON':{'months':[9,10,11],'index':2}, 'DJF':{'months':[12,1,2],'index':3}},overwrite=True):
	"""
	This function reads a state field created by :meth:`temp_anomaly_to_ind` or :meth:`precip_to_index` and finds persistent periods for these statesself. It uses :meth:`optimized_period_identifier`

	Parameters
	----------
		state_file: str
			filepath of a state file. This file needs to have a variable `'state'` with -1 and 1 for the two different states. This file can be created by :meth:`temp_anomaly_to_ind` or :meth:`precip_to_index`
		out_file: str
			filepath of a period file
		seasons: dict, default=`{'MAM':{'months':[3,4,5],'index':0}, 'JJA':{'months':[6,7,8],'index':1}, 'SON':{'months':[9,10,11],'index':2}, 'DJF':{'months':[12,1,2],'index':3}}``
			dictionnary used to cluster detected periods into seasons. If no seasonal analysis is required use `seasons={'year':{'months':range(12),'index':0}}`
		overwrite: bool
			overwrites existing files
	"""

	nc_in=Dataset(state_file,'r')
	# handle time
	time_axis=nc_in.variables['time'][:]
	datevar = num2date(time_axis,units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar)
	month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])
	year=np.array([int(str(date).split("-")[0])	for date in datevar[:]])

	season=month.copy()*np.nan
	for sea in seasons.keys():
		season[np.where((month==seasons[sea]['months'][0]) | (month==seasons[sea]['months'][1]) | (month==seasons[sea]['months'][2]) )[0]]=seasons[sea]['index']

	monthly_index=np.array([mon+yr*12 for mon,yr in zip(month-1,year-np.min(year))])
	mon_year_axis=np.array([yr+mn*0.01 for yr,mn in zip(year,month)])

	state=nc_in.variables['state'][:,:,:]

	#period=state.copy()*np.nan
	period_length=state.copy()*np.nan
	period_state=state.copy()*np.nan
	period_midpoints=state.copy()*np.nan
	period_season=state.copy()*np.nan
	period_monthly_index=state.copy()*np.nan

	period_number=[]
	for y in range(state.shape[1]):
		for x in range(state.shape[2]):
			start_time=time.time()
			try:
				periods=optimized_period_identifier(state[:,y,x].copy())
				identified_periods=np.where(periods!=0)[0]
				per_num=len(identified_periods)
				period_number.append(per_num)

				period_length[0:per_num,y,x]=periods[identified_periods]
				period_state[0:per_num,y,x]=np.sign(periods[identified_periods])
				period_midpoints[0:per_num,y,x]=time_axis[identified_periods]
				period_season[0:per_num,y,x]=season[identified_periods]
				period_monthly_index[0:per_num,y,x]=monthly_index[identified_periods]
			except:
				print('issue at grid ',y,' ',x)

	per_num=max(period_number)

	if overwrite: os.system('rm '+out_file)
	nc_out=Dataset(out_file,'w')
	for dname, the_dim in nc_in.dimensions.iteritems():
		if dname in ['lon','lat']:nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
	nc_out.createDimension('period_id', per_num)

	for v_name, varin in nc_in.variables.iteritems():
		if v_name in ['lon','lat']:
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]

	outVar = nc_out.createVariable('period_length','i2',('period_id','lat','lon',))
	outVar.long_name='period length in days'
	outVar[:] = period_length[0:per_num,:,:]

	outVar = nc_out.createVariable('period_state','i1',('period_id','lat','lon',))
	outVar.long_name='period type'
	outVar[:] = period_state[0:per_num,:,:]

	outVar = nc_out.createVariable('period_midpoints','f',('period_id','lat','lon',))
	outVar.long_name='midpoint of period'
	outVar[:] = period_midpoints[0:per_num,:,:]

	outVar = nc_out.createVariable('period_monthly_index','i2',('period_id','lat','lon',))
	outVar.long_name='monthly index 0 to number of years * 12'
	outVar.first_time_step=str(year[0])+' - '+str(min(month))
	outVar.last_time_step=str(year[-1])+' - '+str(max(month))
	outVar[:] = period_monthly_index[0:per_num,:,:]

	outVar = nc_out.createVariable('period_season','i1',('period_id','lat','lon',))
	outVar.long_name='season in which the midpoint of period is located'
	outVar.description=str(seasons)
	outVar[:] = period_season[0:per_num,:,:]

	nc_out.close()
	nc_in.close()

def temp_anomaly_to_ind(anom_file,out_file,var_name='tas',seasons={'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11],'DJF':[12,1,2]},overwrite=True):
	"""
	Classifies daily temperature anomalies into 'cold' and 'warm' days using the season and grid-cell specific median as threshold

	Parameters
	----------
		anom_file: str
			filepath of a temperature anomalies file. The variable that is read in can be specified with `var_name`.
		out_file: str
			filepath of a state file
		var_name: str
			name of the variable read in `anom_file`
		seasons: dict, default=`{'MAM':{'months':[3,4,5],'index':0}, 'JJA':{'months':[6,7,8],'index':1}, 'SON':{'months':[9,10,11],'index':2}, 'DJF':{'months':[12,1,2],'index':3}}``
			dictionnary used to cluster detected periods into seasons. If no seasonal analysis is required use `seasons={'year':{'months':range(12),'index':0}}`
		overwrite: bool
			overwrites existing files
	"""
	nc_in=Dataset(anom_file,'r')
	time=nc_in.variables['time'][:]
	datevar = num2date(time,units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar)
	month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])

	anom=nc_in.variables[var_name][:,:,:]

	for season in seasons.keys():
		days_in_season=np.where( (month==seasons[season][0]) | (month==seasons[season][1]) | (month==seasons[season][2]) )[0]
		seasonal_median=np.nanmedian(anom[days_in_season,:,:],axis=0)
		anom[days_in_season,:,:]-=seasonal_median

	anom[anom>=0] = 1
	anom[anom<0] = -1

	if overwrite: os.system('rm '+out_file)
	nc_out=Dataset(out_file,'w')
	for dname, the_dim in nc_in.dimensions.iteritems():	nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
	for v_name, varin in nc_in.variables.iteritems():
		if v_name!=var_name:
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]
		else:
			outVar = nc_out.createVariable('state','i1',('time','lat','lon',))
			outVar.description='daily anomalies - seasonal medain of daily anomalies at grid cell level. positive anomalies -> 1 negative anomalies -> -1'
			outVar[:] = anom

	nc_out.close()
	nc_in.close()

def precip_to_index(in_file,out_file,var_name='pr',unit_multiplier=1,threshold=0.5,overwrite=True):
	"""
	Classifies daily precipitation into 'wet' and 'dry' days based on a `threshold`

	Parameters
	----------
		anom_file: str
			filepath of a daily precipitation file. The variable that is read in can be specified with `var_name`.
		out_file: str
			filepath of a state file
		var_name: str
			name of the variable read in `anom_file`
		threshold: float,default=0.5
			threshold used to differentiate between wet and dry days
		unit_multiplier: float,default=1
			factor to multiply daily precipiation with to get mm as units
		overwrite: bool
			overwrites existing files
	"""
	nc_in=Dataset(in_file,'r')
	anom=np.ma.getdata(nc_in.variables[var_name][:,:,:])*unit_multiplier
	mask=np.ma.getmask(nc_in.variables[var_name][:,:,:])
	anom[mask]=np.nan


	anom[anom>=threshold] = 1
	anom[anom<threshold] = -1
	anom[anom**2!=1]=np.nan

	if overwrite: os.system('rm '+out_file)
	nc_out=Dataset(out_file,'w')
	for dname, the_dim in nc_in.dimensions.iteritems():	nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
	for v_name, varin in nc_in.variables.iteritems():
		if v_name!=var_name:
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]
		else:
			outVar = nc_out.createVariable('state','i1',('time','lat','lon',),fill_value=2)
			outVar.description='dry (wet) days are days with precipiation below (above) '+str(threshold)+' and are saved as -1 (1)'
			outVar[:] = anom

	nc_out.close()
	nc_in.close()
