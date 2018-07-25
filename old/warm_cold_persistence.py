import os,sys,glob,time,collections
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import random as random

os.chdir('/Users/peterpfleiderer/Documents/Projects/persistence_py')
try:del sys.modules['persistence_functions'] 
except:pass
from persistence_functions import *
os.chdir('/Users/peterpfleiderer/Documents/Projects')

# test
# temp_anomaly_to_ind('tas_anom_Aday_MIROC5_Plus15-Future_CMIP5-MMM-est1_v2-0_run001.nc','tas_anom_Aday_MIROC5_Plus15-Future_CMIP5-MMM-est1_v2-0_run001_state.nc')

working_path='/global/cscratch1/sd/pepflei/MIROC/'

# # test 
# in_file='/global/cscratch1/sd/pepflei/MIROC/Plus15-Future/tas_Aday_MIROC5_Plus15-Future_CMIP5-MMM-est1_v2-0_run069_land_state.nc'
# get_persistence(in_file,in_file.replace('_state.nc','_period.nc'),overwrite=True)

# adasdasd

scenario_list=[path.split('/')[-1] for path in glob.glob(working_path+'*')]
for scenario in ['Plus20-Future']:
	all_files=glob.glob(working_path+scenario+'/*')
	for file in all_files:
		if len(file.split('_'))==7:
			print file
			start_time=time.time()

			# step 2
			in_file=file
			out_file=in_file.replace('.nc','_land.nc')
			os.system('cdo mul '+in_file+' /global/homes/p/pepflei/masks/landmask_128x256_NA-1.nc '+out_file)

			# step 3
			in_file=out_file
			a=in_file.replace('.nc','_a.nc')
			b=in_file.replace('.nc','_b.nc')
			os.system('cdo trend '+in_file+' '+a+' '+b)
			detrend_1=in_file.replace('.nc','_detrend_1.nc')
			os.system('cdo subtrend '+in_file+' '+a+' '+b+' '+detrend_1)
			os.system('rm '+a+' '+b)

			runmean=in_file.replace('.nc','_runmean.nc')
			os.system('cdo runmean,90 '+detrend_1+' '+runmean)

			detrend_cut=in_file.replace('.nc','_detrend_cut.nc')
			command='cdo delete,timestep='
			for i in range(1,46,1): command+=str(i)+','
			for i in range(1,46,1): command+=str(-i)+','
			os.system(command+' '+detrend_1+' '+detrend_cut)
			out_file=in_file.replace('.nc','_anom.nc')
			os.system('cdo sub '+detrend_cut+' '+runmean+' '+out_file)
			# clean
			os.system('rm '+detrend_1+' '+detrend_cut+' '+runmean+' '+in_file)

			# step 4
			in_file=out_file
			out_file=in_file.replace('_anom.nc','_state.nc')	
			temp_anomaly_to_ind(in_file,out_file,overwrite=True)
			os.system('rm '+in_file)

			# step 5
			in_file=out_file
			try:
				get_persistence(in_file,in_file.replace('_state.nc','_period.nc'),overwrite=True)
				print 'processing time:',time.time()-start_time
			except:
				failed_files=open('output/step_5_period_fails.txt','w')
				failed_files.write(in_file+'\n')
				failed_files.close()
				print '/!\---------------/!\ \n failed for',in_file,'\n/!\---------------/!\ '
















