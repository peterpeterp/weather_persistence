
# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

import os,sys,glob,time,collections,gc
import numpy as np
from netCDF4 import Dataset,num2date
import random as random
import dimarray as da

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
		if np.isfinite(ind[i])==False and np.isfinite(state):
			pers[i-count/2-1]=state*count
			state=np.nan
			count=1
		elif np.isfinite(ind[i]) and np.isfinite(state)==False:
			state=ind[i]
			count=1
		elif ind[i]==state*-1:
			pers[i-count/2-1]=state*count
			count=1
			state*=-1
		elif ind[i]==state:
			count+=1

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

	ind_tmp=ind.copy()
	ind_tmp[ind_tmp!=1]=999
	ind_tmp[ind_tmp==1]=0
	ind_tmp[ind_tmp==999]=1
	su=np.cumsum(ind_tmp)
	counter=collections.Counter(su)

	index=0
	for count,val in zip(counter.values(),counter.keys()):
		index+=count
		if count>=1:
			pers[index-(count-1)/2-1]=1*(count-1)
	# correct start
	if len(ind_tmp)==1:	pers[0]=1
	else:
		if ind_tmp[0]==0 and ind_tmp[1]==1:	pers[0]=1
		if ind_tmp[0]==0 and ind_tmp[1]==0:	pers[np.where(pers<0)[0][0]]+=1

	ind_tmp=ind.copy()
	ind_tmp[ind_tmp!=-1]=999
	ind_tmp[ind_tmp==-1]=0
	ind_tmp[ind_tmp==999]=1
	su=np.cumsum(ind_tmp)
	counter=collections.Counter(su)

	index=0
	for count,val in zip(counter.values(),counter.keys()):
		index+=count
		if count>1:
			pers[index-(count-1)/2-1]=-(count-1)
	# correct start
	if len(ind_tmp)==1:	pers[0]=-1
	else:
		if ind_tmp[0]==0 and ind_tmp[1]==1:	pers[0]=-1
		if ind_tmp[0]==0 and ind_tmp[1]==0:	pers[np.where(pers<0)[0][0]]-=1

	return(pers)


def test_persistence(N):
	ind=np.random.random(N)
	ind[np.where((ind<0.6) & (ind>0.4))[0]]=np.nan
	ind[ind<0.5]=-1
	ind[ind>=0.5]=1
	ind[:4]=[1,np.nan,1,np.nan]
	ind=np.array(ind,'f')
	print(ind[0:100])

	start_time = time.time()
	print(period_identifier(ind)[0:100])
	print("--- basic_and_understandable %s seconds ---" % (time.time() - start_time))

	start_time = time.time()
	print(optimized_period_identifier(ind)[0:100])
	print("--- optimized_period_identifier %s seconds ---" % (time.time() - start_time))

#test_persistence(100)

def get_persistence(state_file,out_file, lat_name='lat', lon_name='lon', seasons={'MAM':{'months':[3,4,5],'index':0}, 'JJA':{'months':[6,7,8],'index':1}, 'SON':{'months':[9,10,11],'index':2}, 'DJF':{'months':[12,1,2],'index':3}},overwrite=True):
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
	if 'calendar' in nc_in.variables['time'].ncattrs():
		datevar = num2date(time_axis,units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar)
	else:
		datevar = num2date(time_axis,units = nc_in.variables['time'].units)

	month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])
	year=np.array([int(str(date).split("-")[0])	for date in datevar[:]])

	season=month.copy()*np.nan
	for sea in seasons.keys():
		season[np.where((month==seasons[sea]['months'][0]) | (month==seasons[sea]['months'][1]) | (month==seasons[sea]['months'][2]) )[0]]=seasons[sea]['index']

	monthly_index=np.array([mon+yr*12 for mon,yr in zip(month-1,year-np.min(year))])
	mon_year_axis=np.array([yr+mn*0.01 for yr,mn in zip(year,month)])

	state=np.ma.getdata(nc_in.variables['state'][:].squeeze())
	mask=np.ma.getmask(nc_in.variables['state'][:].squeeze())
	state[mask]=np.nan
	print(state.shape)
	print(np.nanpercentile(state,range(100)))


	#period=state.copy()*np.nan
	period_length=state.copy()*np.nan
	period_state=state.copy()*np.nan
	period_midpoints=state.copy()*np.nan
	period_season=state.copy()*np.nan
	period_monthly_index=state.copy()*np.nan
	gc.collect()

	period_number=[]
	print('finding periods\n10------50-------100')
	for y,progress in zip(range(state.shape[1]), np.array([['-']+['']*(state.shape[1]/20+1)]*20).flatten()[0:state.shape[1]]):
		sys.stdout.write(progress); sys.stdout.flush()
		for x in range(state.shape[2]):
			start_time=time.time()
			if np.isfinite(np.nanmean(state[:,y,x])) and np.nanmean(state[:,y,x]) not in [-1,1,0]:
				print(np.nanmean(state[:,y,x]))
				periods=optimized_period_identifier(state[:,y,x].copy())
				identified_periods=np.where(periods!=0)[0]
				per_num=len(identified_periods)
				period_number.append(per_num)

				period_length[0:per_num,y,x]=periods[identified_periods]
				period_state[0:per_num,y,x]=np.sign(periods[identified_periods])
				period_midpoints[0:per_num,y,x]=time_axis[identified_periods]
				period_season[0:per_num,y,x]=season[identified_periods]
				period_monthly_index[0:per_num,y,x]=monthly_index[identified_periods]
			# except:
			# 	print('issue at grid ',y,' ',x)
			gc.collect()

	per_num=max(period_number)

	if overwrite: os.system('rm '+out_file)
	nc_out=Dataset(out_file,'w')
	for dname, the_dim in nc_in.dimensions.iteritems():
		if dname in [lon_name,lat_name]:nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
	nc_out.createDimension('period_id', per_num)

	for v_name, varin in nc_in.variables.iteritems():
		if v_name in [lon_name,lat_name]:
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
			outVar[:] = varin[:]

	outVar = nc_out.createVariable('period_length','i2',('period_id',lat_name,lon_name,))
	outVar.long_name='period length in days'
	outVar[:] = period_length[0:per_num,:,:]

	outVar = nc_out.createVariable('period_state','i1',('period_id',lat_name,lon_name,))
	outVar.long_name='period type'
	outVar[:] = period_state[0:per_num,:,:]

	outVar = nc_out.createVariable('period_midpoints','f',('period_id',lat_name,lon_name,))
	outVar.long_name='midpoint of period'
	outVar[:] = period_midpoints[0:per_num,:,:]

	outVar = nc_out.createVariable('period_monthly_index','i2',('period_id',lat_name,lon_name,))
	outVar.long_name='monthly index 0 to number of years * 12'
	outVar.first_time_step=str(year[0])+' - '+str(min(month))
	outVar.last_time_step=str(year[-1])+' - '+str(max(month))
	outVar[:] = period_monthly_index[0:per_num,:,:]

	outVar = nc_out.createVariable('period_season','i1',('period_id',lat_name,lon_name,))
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
	nc=da.read_nc(anom_file)
	datevar=num2date(nc['time'].values,units = nc['time'].units)
	month=np.array([date.month for date in datevar])

	anom=nc[var_name][:,:,:]

	state=nc[var_name].copy()*np.nan

	for season in seasons.keys():
		days_in_season=np.where( (month==seasons[season][0]) | (month==seasons[season][1]) | (month==seasons[season][2]) )[0]
		seasonal_median=np.nanmedian(anom.ix[days_in_season,:,:],axis=0)
		anom.ix[days_in_season,:,:]-=seasonal_median

	state[anom>=0] = 1
	state[anom<0] = -1

	if overwrite: os.system('rm '+out_file)
	state.description='daily anomalies - seasonal medain of daily anomalies at grid cell level. positive anomalies -> 1 negative anomalies -> -1'
	da.Dataset({'state':state}).write_nc(out_file)


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
	nc=da.read_nc(in_file)
	pr=np.ma.getdata(nc[var_name][:,:,:])*unit_multiplier
	mask=np.ma.getmask(nc[var_name][:,:,:])
	pr[mask]=np.nan
	# pr=nc[var_name].values*unit_multiplier

	state=nc[var_name].copy()*np.nan

	state[pr>=threshold] = 1
	state[pr<threshold] = -1
	state[state**2!=1]=np.nan

	#state.values=np.array(state.values,np.int)

	if overwrite: os.system('rm '+out_file)
	state.description='dry (wet) days are days with precipiation below (above) '+str(threshold)+' and are saved as -1 (1)'
	da.Dataset({'state':state}).write_nc(out_file)


def compound_precip_temp_index(tas_state_file,pr_state_file,out_file,overwrite=True):
	"""
	Not documented yet
	"""

	# tas_state=da.read_nc(tas_state_file)['state']
	# pr_state=da.read_nc(pr_state_file)['state']
	# print(np.nanpercentile(tas_state,range(100)))
	# print(np.nanpercentile(pr_state,range(100)))

	nc=da.read_nc(tas_state_file)
	tas_state=np.ma.getdata(nc['state'][:,:,:])
	mask=np.ma.getmask(nc['state'][:,:,:])
	tas_state[mask]=np.nan

	nc=da.read_nc(pr_state_file)
	pr_state=np.ma.getdata(nc['state'][:,:,:])
	mask=np.ma.getmask(nc['state'][:,:,:])
	pr_state[mask]=np.nan

	compound_state_tmp = tas_state.copy()+pr_state.copy()*10
	print(np.nanpercentile(compound_state_tmp,range(100)))
	compound_state_tmp[compound_state_tmp==-9] = 1
	compound_state_tmp[compound_state_tmp==9] = -1
	compound_state_tmp[compound_state_tmp**2!=1]=np.nan


	if overwrite: os.system('rm '+out_file)
	compound_state=da.read_nc(tas_state_file)['state'].copy()
	compound_state.values=compound_state_tmp
	compound_state.description='warm-dry (cold-wet) days are saved as 1 (-1)'
	da.Dataset({'state':compound_state}).write_nc(out_file)



# def optimized_period_identifier(ind):
# 	"""
# 	This function identifies persistent periods using collections. It isn't a straight foreward implementation but runs faster than :meth:`period_identifier`
#
# 	Parameters
# 	----------
# 		ind: np.array
# 			array containing -1 and 1 corresponding to two state
#
# 	Returns
# 	--------
# 		pers: np.array
# 			array of the same length as `ind` containing the length of identified periods. Periods of state -1 are have negative lengths. The period length of a period is placed in the center of the period. All other values are 0.
# 	"""
# 	pers=ind.copy()*0
#
# 	ind[ind<0]=0
#
# 	cuts=list(np.where(np.isfinite(ind)==False)[0])
# 	cuts.append(len(ind))
# 	cut_start=0
#
# 	if len(np.where(ind==0)[0])==len(ind):
# 		return pers
#
# 	for cut_stop in cuts:
# 		if cut_start==cut_stop:
# 			cut_start=cut_stop+1
# 		else:
# 			ind_cut=ind[cut_start:cut_stop]
# 			pers_cut=ind_cut.copy()*0
#
# 			su=np.cumsum(ind_cut)
# 			counter=collections.Counter(su)
#
# 			index=0
# 			for count,val in zip(counter.values(),counter.keys()):
# 				index+=count
# 				if count>1:
# 					pers_cut[index-(count-1)/2-1]=-1*(count-1)
# 			# correct start
# 			if len(ind_cut)==1:	pers_cut[0]=-1
# 			else:
# 				if ind_cut[0]==0 and ind_cut[1]==1:	pers_cut[0]=-1
# 				if ind_cut[0]==0 and ind_cut[1]==0:	pers_cut[np.where(pers_cut<0)[0][0]]-=1
#
# 			ind_cut=-ind_cut+1
# 			su=ind_cut.copy()*0 + np.nan_to_num(ind_cut).cumsum()
# 			counter=collections.Counter(su)
#
# 			index=0
# 			for count,val in zip(counter.values(),counter.keys()):
# 				index+=count
# 				if count>1:
# 					pers_cut[index-(count-1)/2-1]=count-1
# 			# correct start
# 			if len(ind_cut)==1:	pers_cut[0]=1
# 			else:
# 				if ind_cut[0]==0 and ind_cut[1]==1:	pers_cut[0]=1
# 				if ind_cut[0]==0 and ind_cut[1]==0:	pers_cut[np.where(pers_cut>0)[0][0]]+=1
#
# 			pers[cut_start:cut_stop]=pers_cut
# 			cut_start=cut_stop+1
#
# 	return(pers)
