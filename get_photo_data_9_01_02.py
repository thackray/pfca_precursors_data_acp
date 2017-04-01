import gchem
import pylab as pl
import numpy as np
from pytools import pdump
from pytools import dictcsv

JJ = dictcsv('photodata.csv')
print JJ['Lat'][90::2]
Jwinter = np.array(JJ['Jan'][90::2])
Jwinter /= max(Jwinter)
Jsummer = np.array(JJ['Jul'][90::2])
Jsummer /= max(Jsummer)

print Jwinter
print Jsummer
print len(Jsummer)

summerscale = []
for i in range(72):
    summerscale.append(Jsummer)
Jsummer = np.array(summerscale)

winterscale = []
for i in range(72):
    winterscale.append(Jwinter)
Jwinter = np.array(winterscale)

fil = 'summer_photo.bpch'
fil2 = 'full.bpch'
bpch_f = gchem.bpch.open_file(fil)
bpch_f2 = gchem.bpch.open_file(fil2)

for i, db in enumerate(bpch_f.datablocks[:]):
    print i, db, np.mean(db.value[:]), np.max(db.value[:])  # print choices a lgamap                                                                         

OHa = bpch_f.datablocks[12] 
OHb = bpch_f.datablocks[29] 
HO2 = bpch_f2.datablocks[439] #439 #729
NOa = bpch_f.datablocks[15]
NO2a = bpch_f.datablocks[13]
NOb = bpch_f.datablocks[32]
NO2b = bpch_f.datablocks[30]
prop = bpch_f.datablocks[22]
eth = bpch_f.datablocks[24]
J = bpch_f.datablocks[18]


Ta = bpch_f.datablocks[16]
Tb = bpch_f.datablocks[33]
Nair = bpch_f.datablocks[14]
Nair = Nair.value[:,:,:]

print np.shape(OHa)
print np.shape(Nair)

lvlmin, lvlmax = 0,1
NOa = NOa.value[:,:,lvlmax]
NOb = NOb.value[:,:,lvlmax]
NO2a = NO2a.value[:,:,lvlmax]
NO2b = NO2b.value[:,:,lvlmax]
OHa = OHa.value[:,:,lvlmax]
OHb = OHb.value[:,:,lvlmax]
Ta = Ta.value[:,:,lvlmax]
Tb = Tb.value[:,:,lvlmax]
NO = np.where(NOa>NOb,NOa,NOb)*Nair[:,:,lvlmax]
NO2 = np.where(NO2a>NO2b,NO2a,NO2b)*Nair[:,:,lvlmax]
OH = np.where(OHa>OHb,OHa,OHb)
HO2 = HO2.value[:,:,lvlmax]*Nair[:,:,lvlmax]
T = np.where(Ta>Tb,Ta,Tb)
eth = eth.value[:,:,lvlmax]*Nair[:,:,lvlmax]
prop = prop.value[:,:,lvlmax]*Nair[:,:,lvlmax]
J = J.value[:,:,lvlmax]
NNN = Nair[:,:,lvlmax]
print np.mean(NNN)

print np.max(eth), np.max(prop)


k1a = 2.45e-12*np.exp(-1775./T)
k1b = 7.66e-12*np.exp(-1020./T)
k1c = 8.7e-12*np.exp(-615./T)
k3 = 2.8e-12*np.exp(300/T)
k5 = 6e-13*np.exp(725/T)
ch4 = 1.8e-6*Nair[:,:,lvlmax]


hv = 1.e15*np.ones_like(J)*Jsummer
print np.max(hv), np.mean(hv), np.min(hv)

RO2 = (ch4*k1a*OH+eth*k1b*OH+prop*k1c*OH)/(NO*k3+HO2*k5)
print np.max(RO2), np.mean(RO2), np.min(RO2)

dat = {'NO':NO, 'HO2':HO2, 'NO2':NO2, 'OH':OH, 'T':T, 'CH3O2':RO2, 'hv350':hv}
pdump(dat, 'summer_photo4.pkl')
for key in dat:
    print np.shape(dat[key])

