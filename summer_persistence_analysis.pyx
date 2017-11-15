from __future__ import division
from numpy cimport ndarray
import numpy as np
cimport numpy as np
cimport cython

ctypedef np.int32_t np_int_t
ctypedef np.float64_t np_float_t

def summer_period_analysis(np.ndarray[np_int_t, ndim=3] per_len, np.ndarray[np_int_t, ndim=3] per_mid, np.ndarray[np_int_t, ndim=3] summer_events, np.ndarray[np_int_t, ndim=2] mask, np.ndarray[np_float_t, ndim=3] tas, np.ndarray[np_int_t, ndim=1] year_var, int Ni, int Nx, int Ny):
  cdef np.ndarray[np_int_t, ndim=3] hot_shift = np.zeros([Ni,Ny,Nx], dtype=np.int32)
  cdef np.ndarray[np_float_t, ndim=3] hot_temp = np.zeros([Ni,Ny,Nx], dtype=np.float64)
  cdef np.ndarray[np_int_t, ndim=3] yearly_X = np.zeros([Ni,Ny,Nx], dtype=np.int32)
  cdef np.ndarray[np_float_t, ndim=3] cum_heat = np.zeros([Ni,Ny,Nx], dtype=np.float64)
  cdef int start_year = year_var[0]
  cdef int end_year = year_var[-1]
  cdef np.ndarray[np_float_t, ndim=3] tasX = np.zeros([Ni,Ny,Nx], dtype=np.float64)
  cdef int Nevents = summer_events.shape[0]
  for i in range(Nevents):
    for y in range(Ny):
      if mask[y,x]==1:
        for yr in range(start_year,end_year):
          for i in range(Ni):
            if per_mid[i,y,x]>0:
              if per_sea[i,y,x]==1 and year_var[per_mid[i,y,x]]==yr:
                tasX[i,y,x]=max(tas[np.where(year_var==yr)[0],y,x])
                days=range(per_mid[i,y,x]-int(abs(per_len[i,y,x])/2.),per_mid[i,y,x]+int(round(abs(per_len[i,y,x])/2.)))
                hot_temp[i,y,x]=np.max(tas[days,y,x])
                cum_heat[i,y,x]=np.sum(tas[days,y,x])
                hot_shift[i,y,x]=days[np.argmax(tas[days,y,x])]-per_mid[i,y,x]
                if hot_temp[i,y,x]==tasX[i,y,x]:
                  yearly_X[i,y,x]=1
  return  cum_heat,hot_shift,hot_temp,yearly_X,tasX

def summer_period_analysis__(np.ndarray[np_int_t, ndim=3] ll, np.ndarray[np_int_t, ndim=3] mm, np.ndarray[np_int_t, ndim=3] sea, np.ndarray[np_int_t, ndim=2] mask, np.ndarray[np_float_t, ndim=3] tt, int Ni, int Nx, int Ny):
  cdef np.ndarray[np_int_t, ndim=3] ho = np.zeros([Ni,Ny,Nx], dtype=np.int32)
  cdef np.ndarray[np.float_t, ndim=3] cu = np.zeros([Ni,Ny,Nx], dtype=np.float64)
  for x in range(Nx):
    for y in range(Ny):
      for i in range(Ni):
        if sea[i,y,x]==1 and mask[y,x]==1:
          days=range(mm[i,y,x]-int(abs(ll[i,y,x])/2.),mm[i,y,x]+int(round(abs(ll[i,y,x])/2.)))
          cu[i,y,x]=np.nansum(tt[days,y,x])
          ho[i,y,x]=days[np.argmax(tt[days,y,x])]-mm[i,y,x]
  return  cu,ho
