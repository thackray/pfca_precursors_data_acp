import os
import sys
from pytools import pdump, pload
from datetime import datetime, timedelta
from boxmodel import BoxModel
from hysplit import getGCbox_fromfile
import numpy as np


class ThisModel(BoxModel):
    def set_tout(self,tout):
        self.chembox.meta.set_tout(tout)
        print(self.chembox.meta.get_constant_species())
        return
    def set_HO2(self, val):
        self.chembox.concs['HO2'] = val
    def set_NO(self, val):
        self.chembox.concs['NO'] = val
    def set_XO2(self, val):
        self.chembox.concs['XO2'] = val
    def set_NO2(self, val):
        self.chembox.concs['NO2'] = val
    def set_OH(self, val):
        self.chembox.concs['OH'] = val
    def set_T(self, val):
        self.chembox.inputs.conditions['T'] = val
    def set_conc(self, spc, val):
        self.chembox.concs[spc] = val

rundir = ''
hyI,hyJ = getGCbox_fromfile('hysplit_traj.txt')
times = range(len(hyI))
dic = pload(os.path.join(rundir,'summer_photo.pkl'),encoding="bytes")
ii = 0
res_grid = []
B = ThisModel(os.path.join(rundir, 'config.yaml'),
              os.path.join(rundir, 'inputs.yaml'))


diccy = dic.copy()
for key in diccy:
    diccy[key] = []
for t,T in enumerate(times):
    II = hyI[t]
    JJ = hyJ[t]
    for key in diccy:
        diccy[key].append(dic[key][II, JJ])

for key in diccy:
    diccy[key] = np.interp(range(len(times)),range(0,len(times),3),diccy[key][::3])

outout = {}

for t,T in enumerate(times):
    
    tout = [datetime(1987,6,7),datetime(1987,6,7,1)]
    B.set_tout(tout)
    for key in diccy:
        if key == 'T':
            print(diccy['T'][t])
            B.set_T(diccy['T'][t])
        else:
            print(key, diccy[key][t])
            B.set_conc(key,diccy[key][t])
    output = B.run()
    if t == 0:
        for key in output['concs']:
            outout[key] = []
    for key in output['concs']:
        outout[key].append(output['concs'][key][-1])
    print(output['concs']['C8F17CH2C(O)H'][-1])
    print(output['concs']['C8F17C(O)OH'][-1])
    print(output['concs']['C7F15C(O)OH'][-1])

pdump(outout,'parcelout.pkl')
