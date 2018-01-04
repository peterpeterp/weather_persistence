from __future__ import division
from numpy cimport ndarray
import numpy as np
cimport numpy as np
cimport cython

ctypedef np.int32_t np_int_t
ctypedef np.float64_t np_float_t

def summer_period_analysis(np.ndarray[np_int_t, ndim=3] per_len, np.ndarray[np_int_t, ndim=3] per_mid, np.ndarray[np_int_t, ndim=3] seas, np.ndarray[np_int_t, ndim=3] state, np.ndarray[np_float_t, ndim=3] tas, np.ndarray[np_float_t, ndim=2] x90_thresh, np.ndarray[np_int_t, ndim=1] year_var, int Ny, int Nx, int Ni):

  cdef np.ndarray[np_int_t, ndim=2] events_in_loc = np.zeros([Ny,Nx], dtype=np.int32)
  for x in range(Nx):
    for y in range(Ny):
      for i in range(Ni):
        if seas[i,y,x]==1 and state[i,y,x]==1 and per_len[i,y,x]>=x90_thresh[y,x]:
          events_in_loc[y,x]+=1

  cdef int Ni_new = np.max(events_in_loc)

  cdef np.ndarray[np_int_t, ndim=3] hot_shift = np.zeros([Ni_new,Ny,Nx], dtype=np.int32)
  cdef np.ndarray[np_float_t, ndim=3] hot_temp = np.zeros([Ni_new,Ny,Nx], dtype=np.float64)
  cdef np.ndarray[np_float_t, ndim=3] cum_heat = np.zeros([Ni_new,Ny,Nx], dtype=np.float64)
  cdef np.ndarray[np_float_t, ndim=3] tasX = np.zeros([Ni_new,Ny,Nx], dtype=np.float64)
  cdef np.ndarray[np_int_t, ndim=3] original_period_id = np.zeros([Ni_new,Ny,Nx], dtype=np.int32)

  print Ni_new

  cdef int count=0
  for x in range(Nx):
    for y in range(Ny):
      count=0
      for i in range(Ni):
        if seas[i,y,x]==1 and state[i,y,x]==1 and per_len[i,y,x]>=x90_thresh[y,x]:
          tasX[count,y,x]=np.max(tas[np.where(year_var==year_var[per_mid[i,y,x]])[0],y,x])
          low=int(abs(per_len[i,y,x])/2.)
          high=int(round(abs(per_len[i,y,x])/2.))
          shift=low%2
          days=range(per_mid[i,y,x]-low+shift,per_mid[i,y,x]+high+shift)
          #days=range(per_mid[i,y,x]-int(abs(per_len[i,y,x])/2.),per_mid[i,y,x]+int(round(abs(per_len[i,y,x])/2.)))
          original_period_id[count,y,x]=i
          hot_temp[count,y,x]=np.max(tas[days,y,x])
          cum_heat[count,y,x]=np.sum(tas[days,y,x])
          hot_shift[count,y,x]=days[np.argmax(tas[days,y,x])]-per_mid[i,y,x]
          count=count+1

  return  cum_heat,hot_shift,hot_temp,tasX,Ni_new,original_period_id
