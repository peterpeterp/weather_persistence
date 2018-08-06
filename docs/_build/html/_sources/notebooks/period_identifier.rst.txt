
Methods used to identify periods
--------------------------------

.. code:: ipython2

    import os,sys,glob,time,collections
    import numpy as np
    from netCDF4 import Dataset,num2date
    import random as random
    import dimarray as da
    sys.path.append(os.path.abspath('../../'))
    from persistence_functions import *

.. code:: ipython2

    ind=np.random.random(100)
    ind[ind<0.5]=-1
    ind[ind>=0.5]=1#
    ind[-2]=np.nan
    ind[5]=np.nan
    ind=np.array(ind,'f')
    print 'state index:',ind[0:100] 


.. parsed-literal::

    state index: [-1. -1.  1. -1. -1. nan  1.  1. -1.  1.  1. -1.  1.  1. -1.  1.  1.  1.
      1. -1.  1.  1.  1. -1.  1.  1.  1.  1. -1. -1.  1. -1. -1.  1.  1. -1.
     -1.  1.  1.  1.  1.  1. -1.  1. -1. -1. -1.  1.  1. -1. -1. -1. -1. -1.
      1.  1. -1. -1. -1. -1.  1. -1.  1. -1.  1. -1.  1.  1.  1.  1.  1. -1.
      1. -1. -1. -1.  1.  1.  1. -1. -1.  1. -1. -1.  1.  1.  1. -1.  1. -1.
      1.  1.  1.  1. -1. -1. -1. -1. nan  1.]


There are two functions identifying persistent periods in state array:
:meth:``period_identifier`` and :meth:``optimized_period_identifier``.
Both generate the same result. :meth:``period_identifier`` is a straight
foreward implementation, :meth:``optimized_period_identifier`` is
anoptimized implementation using ``collections`` and thus running
faster.

.. code:: ipython2

    print 'period_identifier():',period_identifier(ind)[0:100]
    print 'optimized_period_identifier():',optimized_period_identifier(ind)[0:100]


.. parsed-literal::

    period_identifier(): [-2. -0.  1. -2. -0. nan  2.  0. -1.  2.  0. -1.  2.  0. -1.  0.  4.  0.
      0. -1.  0.  3.  0. -1.  0.  4.  0.  0. -2. -0.  1. -2. -0.  2.  0. -2.
     -0.  0.  0.  5.  0.  0. -1.  1. -0. -3. -0.  2.  0. -0. -0. -5. -0. -0.
      2.  0. -0. -4. -0. -0.  1. -1.  1. -1.  1. -1.  0.  0.  5.  0.  0. -1.
      1. -0. -3. -0.  0.  3.  0. -2. -0.  1. -2. -0.  0.  3.  0. -1.  1. -1.
      0.  4.  0.  0. -0. -4. -0. -0. nan  1.]
    optimized_period_identifier(): [ 0. -2.  1. -2.  0. nan  0.  2. -1.  2.  0. -1.  2.  0. -1.  0.  4.  0.
      0. -1.  0.  3.  0. -1.  0.  4.  0.  0. -2.  0.  1. -2.  0.  2.  0. -2.
      0.  0.  0.  5.  0.  0. -1.  1.  0. -3.  0.  2.  0.  0.  0. -5.  0.  0.
      2.  0.  0. -4.  0.  0.  1. -1.  1. -1.  1. -1.  0.  0.  5.  0.  0. -1.
      1.  0. -3.  0.  0.  3.  0. -2.  0.  1. -2.  0.  0.  3.  0. -1.  1. -1.
      0.  4.  0.  0.  0. -4.  0.  0. nan  1.]


.. parsed-literal::

    /Users/peterpfleiderer/Projects/Persistence/weather_persistence/persistence_functions.py:65: RuntimeWarning: invalid value encountered in less
      ind[ind<0]=0

