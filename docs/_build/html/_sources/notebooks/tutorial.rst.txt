
Tutorial
--------

.. code:: ipython2

    import sys,os
    import numpy as np
    from netCDF4 import Dataset,num2date
    import random as random
    import dimarray as da
    sys.path.append(os.path.abspath('../../'))
    from persistence_functions import *

**Load detrended temperature anomalies**

.. code:: ipython2

    tas=da.read_nc('sample_tas.nc')['tas']
    print 'percentiles',np.nanpercentile(tas.ix[:,0,0].values,[10,50,90])


.. parsed-literal::

    percentiles [-1.65030825  0.01166975  1.68142478]


**convert tas anomalies to a state array**

.. code:: ipython2

    temp_anomaly_to_ind('sample_tas.nc','sample_tas_state.nc',var_name='tas',overwrite=True)
    state=da.read_nc('sample_tas_state.nc')['state']
    print 'temperature anomaly:',tas.ix[:20,0,0].values
    print 'state:',state.ix[:20,0,0].values


.. parsed-literal::

    temperature anomaly: [-1.2905025  -2.1703892  -2.124154   -1.8810792  -1.5038338  -2.1755633
     -1.7776151  -0.68114567  0.17394018  0.6489425  -0.35136604  0.46413946
      0.979579    2.148981    0.6392119  -0.9649811  -0.7494569  -0.5388527
     -0.22809458 -1.2265306 ]
    state: [-1. -1. -1. -1. -1. -1. -1. -1.  1.  1.  1.  1.  1.  1.  1. -1. -1. -1.
      1. -1.]


**compute persistence**

.. code:: ipython2

    get_persistence('sample_tas_state.nc','sample_tas_period.nc',overwrite=True)
    period=da.read_nc('sample_tas_period.nc')
    print 'state:',state.ix[:26,0,0].values
    print 'period length:',period['period_length'].ix[:5,0,0].values
    print 'period midpont:',period['period_midpoints'].ix[:5,0,0].values


.. parsed-literal::

    state: [-1. -1. -1. -1. -1. -1. -1. -1.  1.  1.  1.  1.  1.  1.  1. -1. -1. -1.
      1. -1. -1. -1. -1. -1. -1.  1.]
    period length: [-8.  7. -3.  1. -6.]
    period midpont: [46435.5 46442.5 46447.5 46449.5 46452.5]

