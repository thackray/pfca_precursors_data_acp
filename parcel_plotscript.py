from pytools import pload
import pylab as pl

PFNAout = [0.]
PFOAout = [0.]
PREout = [1.]
#INTout = [0.]
SHORTout = [0.]

outout = pload('parcelout.pkl')
PREout += outout['C8F17CH2C(O)H']
PFNAout += outout['C8F17C(O)OH']
PFOAout += outout['C7F15C(O)OH']
SHORTout += outout['C8F17O']
sumarray = pl.zeros(len(PFNAout))

notinter = ['C8F17O','C7F15C(O)OH','C8F17C(O)OH','C8F17CH2C(O)H','H2O(l)','CH3O2','HO2','NO','NO2','Cl','OH','hv350']
for key in outout:
    if key not in notinter:
        print(key)
        sumarray += pl.array([0.]+outout[key])
INTout = sumarray

pl.figure()
pl.subplot(2,1,1)
pl.plot(PREout,lw=3,label='Precursor')
pl.plot(SHORTout, lw=3,label='Short-chain compounds')
pl.plot(INTout, lw=3, label='Intermediates')
pl.legend(loc='upper center')
pl.xlabel('Days',fontsize=20)
x = range(0,len(PREout),24)
pl.xticks(x,[int(i/24) for i in x])
pl.ylabel('Relative concentration',fontsize=20)
pl.xlim(0,len(PREout))

pl.subplot(2,1,2)
pl.plot(PFNAout,lw=3,label='PFNA')
pl.plot(PFOAout, lw=3,label='PFOA')
pl.legend(loc='upper center')
pl.xlabel('Days',fontsize=20)
pl.ylabel('Relative concentration',fontsize=20)
pl.xticks(x,[int(i/24) for i in x])
pl.xlim(0,len(PREout))
pl.ylim(0.,0.015)

pl.show()