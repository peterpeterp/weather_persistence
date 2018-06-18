import os,sys,glob,time,collections
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import random as random
import dimarray as da


def period_identifier(ind):
	pers=ind.copy()*0

	state,count=ind[0],1
	for i in range(1,len(ind)):
		if ind[i]==state:
			count+=1
		if ind[i]!=state:
			pers[i-count/2-1]=state*count
			count=0
			if ind[i]!=99:
				state*=-1
				count=1

	# still an issue with last spell??
	if state==1:	pers[i-count/2]=state*count
	if state==-1:	pers[i-count/2]=state*count

	return(pers)


def optimized_period_identifier(ind):
	# not straight forward but faster
	# works with nans
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
	# ind[-5:]=-1
	ind[-2]=np.nan
	ind[ind<0.5]=-1
	ind[ind>=0.5]=1
	ind=np.array(ind,'i')
	print(ind[0:100])

	start_time = time.time()
	print(period_identifier(ind)[0:100])
	print("--- basic_and_understandable %s seconds ---" % (time.time() - start_time))

	start_time = time.time()
	print(optimized_period_identifier(ind)[0:100])
	print("--- optimized_period_identifier %s seconds ---" % (time.time() - start_time))

#test_persistence(100)

def get_persistence(state_file,out_file,seasons={'MAM':{'months':[3,4,5],'index':0}, 'JJA':{'months':[6,7,8],'index':1}, 'SON':{'months':[9,10,11],'index':2}, 'DJF':{'months':[12,1,2],'index':3}},overwrite=True,eke_file=None,spi_file=None):

	nc_in=Dataset(state_file,'r')
	# handle time
	time_axis=nc_in.variables['time'][:]
	datevar = num2date(time_axis,units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar)
	month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])
	year=np.array([int(str(date).split("-")[0])	for date in datevar[:]])

	season=month.copy()*np.nan
	for sea in seasons.keys():
		season[np.where((month==seasons[sea]['months'][0]) | (month==seasons[sea]['months'][1]) | (month==seasons[sea]['months'][2]) )[0]]=seasons[sea]['index']

	#monthly_index=np.array([mon+yr*12 for mon,yr in zip(month-1,year-np.min(year))])
	mon_year_axis=np.array([yr+mn*0.01 for yr,mn in zip(year,month)])

	state=nc_in.variables['state'][:,:,:]

	#period=state.copy()*np.nan
	period_length=state.copy()*np.nan
	period_state=state.copy()*np.nan
	period_midpoints=state.copy()*np.nan
	period_season=state.copy()*np.nan

	if eke_file is not None:
		eke=da.read_nc(eke_file)
		datevar = num2date(eke['time'].values,units = eke['time'].attrs['units'],calendar = eke['time'].attrs['calendar'])
		month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])
		year=np.array([int(str(date).split("-")[0])	for date in datevar[:]])
		EKE=da.DimArray(axes=[np.unique(mon_year_axis),eke['lat'].values,eke['lon'].values],dims=['time','lat','lon'])
		for my_i in np.unique(mon_year_axis):
			yr=int(my_i)
			mth=int((my_i-yr)*100)+1
			index=np.where((month==mth) & (year==yr))[0]
			if len(index)==1:
				EKE[my_i,:,:]=eke['EKE'].values.squeeze()[index,:,:]
		period_eke=state.copy()*np.nan
		print(EKE.shape)

	if spi_file is not None:
		spi=da.read_nc(spi_file)
		if 'calendar' in spi['time'].attrs.keys():
			calendar=spi['time'].attrs['calendar']
		else:
			calendar='365_day'
		datevar = num2date(spi['time'].values,units = spi['time'].attrs['units'],calendar = calendar)
		month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])
		year=np.array([int(str(date).split("-")[0])	for date in datevar[:]])
		SPI=da.DimArray(axes=[np.unique(mon_year_axis),spi['lat'].values,spi['lon'].values],dims=['time','lat','lon'])
		for my_i in np.unique(mon_year_axis):
			yr=int(my_i)
			mth=int((my_i-yr)*100)+1
			index=np.where((month==mth) & (year==yr))[0]
			if len(index)==1:
				SPI[my_i,:,:]=spi['SPI'].values.squeeze()[index,:,:]
		period_spi=state.copy()*np.nan
		print(SPI.shape)

	period_number=[]
	for y in range(state.shape[1]):
		print(y)
		for x in range(state.shape[2]):
			periods=optimized_period_identifier(state[:,y,x].copy())

			identified_periods=np.where(periods!=0)[0]
			per_num=len(identified_periods)
			period_number.append(per_num)

			period_length[0:per_num,y,x]=periods[identified_periods]
			period_state[0:per_num,y,x]=np.sign(periods[identified_periods])
			period_midpoints[0:per_num,y,x]=time_axis[identified_periods]
			period_season[0:per_num,y,x]=season[identified_periods]
			if eke_file is not None:
				period_eke[0:per_num,y,x]=EKE[mon_year_axis[identified_periods],:,:].ix[:,y,x]
			if spi_file is not None:
				period_spi[0:per_num,y,x]=SPI[mon_year_axis[identified_periods],:,:].ix[:,y,x]

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

	if EKE is not None:
		outVar = nc_out.createVariable('period_eke','f',('period_id','lat','lon',))
		outVar.long_name='monthly EKE of period midpoint'
		outVar[:] = period_eke[0:per_num,:,:]

	if SPI is not None:
		outVar = nc_out.createVariable('period_spi','f',('period_id','lat','lon',))
		outVar.long_name='monthly SPI of period midpoint'
		outVar[:] = period_spi[0:per_num,:,:]

	outVar = nc_out.createVariable('period_season','i1',('period_id','lat','lon',))
	outVar.long_name='season in which the midpoint of period is located'
	outVar.description=str(seasons)
	outVar[:] = period_season[0:per_num,:,:]

	nc_out.close()
	nc_in.close()


def temp_anomaly_to_ind(anom_file,out_file,var_name='tas',seasons={'MAM':[3,4,5],'JJA':[6,7,8],'SON':[9,10,11],'DJF':[12,1,2]},overwrite=True):
	nc_in=Dataset(anom_file,'r')
	time=nc_in.variables['time'][:]
	datevar = num2date(time,units = nc_in.variables['time'].units,calendar = nc_in.variables['time'].calendar)
	month=np.array([int(str(date).split("-")[1])	for date in datevar[:]])

	anom=nc_in.variables[var_name][:,:,:]


	for season in seasons.keys():
		days_in_seayson=np.where( (month==seasons[season][0]) | (month==seasons[season][1]) | (month==seasons[season][2]) )[0]
		seasonal_median=np.nanmedian(anom[days_in_seayson,:,:],axis=0)
		anom[days_in_seayson,:,:]-=seasonal_median

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
