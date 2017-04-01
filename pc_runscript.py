import numpy as np
import sys
import os
import pylab as pl
from datetime import datetime, timedelta

from boxmodel import BoxModel
import mypcm as pcm
from mypc import PCE
from pc_analysis import show_param_fracs, top_coeffs, get_vars, get_skew
from calculate_rates_yields import calculate_timescale as t_y
from pytools import pload, pdump, savefig

class ThisModel(BoxModel):
    def set_tout(self,tout):
        self.chembox.meta.set_tout(tout)
        return

n = 39; o = 2
Pt = PCE(nparams=n, order=o)
Py = PCE(nparams=n, order=o)

inputs = pcm.get_collocation_points(Pt)

rundir = sys.argv[1]

pert_reactions = range(1,40)
pert_reactions = [str(x) for x in pert_reactions]

pfoa, pfna = 'C7F15C(O)OH','C8F17C(O)OH'
schain = 'C8F17O'
spc_interest = [pfoa, pfna, schain]

print np.shape(inputs)
redo_list = []

for ri,inp in enumerate(inputs):
    if (not os.path.exists(os.path.join(rundir,'ensemble/output%s.pkl'%ri))) or\
       (ri in redo_list):
        M = ThisModel(os.path.join(rundir,'config.yaml'),
                      os.path.join(rundir,'inputs.yaml'))
        pert = {}
        for i,p in enumerate(inp):
            pert[pert_reactions[i]] = p
        print ri
#        print inp
#        print pert
        M.rate_perturbations(pert)
        tout = [datetime(1987,6,7)]
        for x in range(1,365):
            tout.append(tout[-1]+timedelta(hours=24))
        M.set_tout(tout)
        output = M.run_to_eq(['C7F15C(O)OH','C8F17O'])
        pdump(output, os.path.join(rundir,'ensemble/output%s.pkl'%ri))
oupt = {}
oupy = {}
for spc in spc_interest:
    oupt[spc] = []
    oupy[spc] = []
probi=[]
for ri,inp in enumerate(inputs):
    o = pload(os.path.join(rundir,'ensemble/output%s.pkl'%ri))
    for spc in spc_interest:
        t1, t2, y = t_y(o['t'], o['concs'][spc])
        if spc is pfna:
            blah, blah, y1 = t_y(o['t'], o['concs']['C8F17CH2OH'])
            blah, blah, y2 = t_y(o['t'], o['concs']['C8F17CHOHOH'])
            y = y
        if (spc == 'C7F15C(O)OH') and (y > 1.):
            print spc, y, ri
            probi.append(ri)
        if (spc == 'C8F17C(O)OH') and (y > 1.):
            print spc, y, ri
            if ri not in probi:
                probi.append(ri)
        if (spc == 'C8F17O') and (y > 1.):
            print spc, y, ri
            if ri not in probi:
                probi.append(ri)
         
        oupt[spc].append(t2.total_seconds()/(24*3600))
        oupy[spc].append(np.log10(y))
print probi

for spc in spc_interest:
    pl.figure()
    pl.hist(oupy[spc])
    pl.title(spc)
    Pt.set_coeffs(pcm.get_coeffs(inputs, oupt[spc], Pt.poly))
    Py.set_coeffs(pcm.get_coeffs(inputs, oupy[spc], Py.poly))
    pdump(Pt, os.path.join(rundir,'timing_poly_%s.pkl'%spc))
    pdump(Py, os.path.join(rundir,'yield_poly_%s.pkl'%spc))

pl.show()
