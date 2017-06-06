import os,sys,glob,time,collections
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import random as random

os.chdir('/Users/peterpfleiderer/Documents/Projects/persistence_py')
try:del sys.modules['persistence_functions'] 
except:pass
from persistence_functions import *
os.chdir('/Users/peterpfleiderer/Documents/Projects/persistence_py')

in_file='test/pr_bced_1960_1999_hadgem2-es_rcp2p6_2011-2020.nc4'

out_file=in_file.replace('.nc4','_dry_wet_state.nc4')	
precip_to_index(in_file,out_file,unit_multiplier=86400,overwrite=True)

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
















