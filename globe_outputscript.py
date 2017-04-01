import os
import sys
from pytools import pdump
from datetime import datetime, timedelta
from boxmodel import BoxModel
import numpy as np
from pytools import pload,pdump
from calculate_rates_yields import calculate_timescale as t_y

pfoa, pfna = 'C7F15C(O)OH','C8F17C(O)OH'
schain = 'C8F17O'
spc_interest = [pfoa, pfna, schain]

rundir = sys.argv[1]
redo = []
ii = 0
res_grid = []
oupt = {'coord':[]}
oupy = {'coord':[]}

for spc in spc_interest:
    oupt[spc] = []
    oupy[spc] = []

mech = 'ensemble_summer_highres'
ii = 0
for II in range(0,72,1):
    for JJ in range(23,46,1):
        print II,JJ
        #if JJ>22:
        if 1:
            o = pload(os.path.join(rundir,'%s/output%s.pkl'%(mech,ii)))
            for spc in spc_interest:
                t1, t2, y = t_y(o['t'], o['concs'][spc])
                oupy[spc].append(y)
                oupt[spc].append(t2)
                oupy['coord'].append((II,JJ))
                oupt['coord'].append((II,JJ))
        ii+=1
#print np.shape(oupy[pfna])
print np.shape(oupt[pfoa])
pdump(oupy,os.path.join(rundir,'%s_yield_output.pkl'%mech))
pdump(oupt,os.path.join(rundir,'%s_timing_output.pkl'%mech))

