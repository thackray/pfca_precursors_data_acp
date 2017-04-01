import os
import sys
from pytools import pdump, pload
from datetime import datetime, timedelta
from boxmodel import BoxModel
import numpy as np


class ThisModel(BoxModel):
    def set_tout(self,tout):
        self.chembox.meta.set_tout(tout)
        print self.chembox.meta.get_constant_species()
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

rundir = sys.argv[1]

dic = pload(os.path.join(rundir,'summer_photo.pkl'))
print 'Running model for '+rundir
redo = []
redo = ['all']
ii = 0
res_grid = []
for II in range(0,72,1):
    for JJ in range(23,46,1):
        print II,JJ
        B = ThisModel(os.path.join(rundir,'config.yaml'), 
                      os.path.join(rundir,'inputs.yaml'))
        tout = [datetime(1987,6,7)]
        for i in range(1,365):
            tout.append(tout[-1]+timedelta(hours=24))
        B.set_tout(tout)
        for key in dic:
            if key == 'T':
                print dic['T'][II,JJ]
                B.set_T(dic['T'][II,JJ])
            else:
                print key, dic[key][II,JJ]
                B.set_conc(key,dic[key][II,JJ])
        if (ii in redo) or ('all' in redo):
            output = B.run_to_eq(['C7F15C(O)OH','C8F17O'])
            pdump(output, os.path.join(rundir,'ensemble_summer_highres/output%s.pkl'%ii))
            print output['concs']['C8F17C(O)OH'][-1]
            print output['concs']['C7F15C(O)OH'][-1]
        ii += 1

