import pylab as pl
from pytools import pload, pdump
import os
import sys
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable as mal

rundir = sys.argv[1]
tag = ''
tag = '_winter'
dic = pload(os.path.join(rundir,'winter_photo.pkl'))

mech = 'ensemble%s'%tag

Y = pload(os.path.join(rundir,'%s_yield_output.pkl'%(mech)))
T = pload(os.path.join(rundir,'%s_timing_output.pkl'%(mech)))

pfoaY = np.array(Y['C7F15C(O)OH'])
pfoaT = np.array([t.days for t in T['C7F15C(O)OH']])
pfnaY = np.array(Y['C8F17C(O)OH'])
pfnaT = np.array([t.days for t in T['C8F17C(O)OH']])
scY = np.array(Y['C8F17O'])
scT = np.array([t.days for t in T['C8F17O']])
summ = pfoaY + pfnaY + scY
coords = np.array(Y['coord'])
print coords
#Z = np.array([c[2] for c in coords[::3]])
Y = np.array([c[1] for c in coords[::3]])
X = np.array([c[0] for c in coords[::3]])
#CVAR = 'OH'
cs = []
for ci,CVAR in enumerate(['OH','HO2','NO','CH3O2','hv350']):
    cs.append(np.array([dic[CVAR][i,j] for i,j in coords[::3]]))
#colorvar = X
colorvar = np.array([dic['NO'][i,j] for i,j in 
                     coords[::3]])
for c in cs:
    print np.mean(c)
print np.shape(X), np.shape(Y), np.shape(pfoaY), np.shape(cs[0])
stack = np.vstack((X,Y,pfnaY,pfnaT,pfoaY,pfoaT,scY,scT,np.log10(cs[0]),
                   np.log10(cs[1]),np.log10(cs[2]),np.log10(cs[3]),
                   np.log10(cs[4])))
pdump(stack, os.path.join(rundir,'%s_data_stack.pkl'%mech))
pdump(['x','y','Y9','t9','Y8','t8','Ys','ts','OH','HO2','NO','RO2',
       'hv350'],os.path.join(rundir,'%s_data_names.pkl'%mech))

ptm = [4,5,6,7] # plots to make

alp = 0.5

fignum = 1
if fignum in ptm:
    pl.figure(fignum)
    pl.subplot(2,2,1)
    pl.scatter(pfoaT,pfoaY, alpha=alp)
    pl.ylabel('PFOA Yield')
    pl.subplot(2,2,2)
    pl.scatter(pfnaT,pfnaY, alpha=alp)
    pl.ylabel('PFNA Yield')
    pl.subplot(2,2,3)
    pl.scatter(scT,scY, alpha=alp)
    pl.ylabel('SC Yield')
    pl.subplot(2,2,4)
    pl.scatter(scT,summ, alpha=alp)

fignum = 2
if fignum in ptm:
    pl.figure(fignum,figsize=(14,8))
    pl.subplot(1,2,1)
    pl.scatter(pfoaT,pfoaY,c=np.log10(colorvar), alpha=alp)
    pl.ylabel('PFOA Yield', fontsize=20)
    pl.xlabel('Formation time (days)', fontsize=20)
    pl.xlim(0,700)
    pl.ylim(0,0.9)
    #pl.colorbar()

    pl.subplot(1,2,2)
    pl.scatter(pfnaT,pfnaY,c=np.log10(colorvar), alpha=alp)
    pl.ylabel('PFNA Yield', fontsize=20)
    pl.xlabel('Formation time (days)', fontsize=20)
    pl.xlim(0,700)
    pl.ylim(0,0.9)
    pl.colorbar()

fignum = 3
if fignum in ptm:
    pl.figure(fignum)
    pl.scatter(colorvar,pfoaY,alpha=0.3)
    pl.semilogx()
    pl.ylim(0.,.4)


left =[0.,0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]#np.linspace(0.,0.9,21)
barleft = range(11)
fignum = 4
if fignum in ptm:
    pl.figure(999, figsize=(13,9))
    n,bins,p=pl.hist(pfoaY, left, color='gray')
    pl.figure(fignum, figsize=(13,9))
    pl.bar(barleft,n, color='gray', width=1,)
    pl.xlabel('PFOA Yield', fontsize=30)
    pl.ylabel('Count', fontsize=30)
    pl.xlim(0,1)
    pl.ylim(0,1000)
    pl.yticks(fontsize=25)
    pl.xticks(barleft,[str(x) for x in left[:-1]],fontsize=25)
    pl.savefig('pfoa_hist%s.png'%tag)

fignum = 5
if fignum in ptm:
    pl.figure(999, figsize=(13,9))
    n,bins,p=pl.hist(pfnaY, left, color='gray')
    pl.figure(fignum, figsize=(13,9))
    #    pl.hist(pfnaY, left, color='gray')
    pl.bar(barleft,n,color='gray',width=1)
    pl.xlabel('PFNA Yield', fontsize=30)
    pl.ylabel('Count', fontsize=30)
    pl.xlim(0,1)
    pl.ylim(0,1000)
    pl.xticks(barleft,[str(x) for x in left[:-1]], fontsize=25)
    pl.yticks(fontsize=25)
    pl.savefig('pfna_hist%s.png'%tag)

X = list(set(X))
Y = list(set(Y))
X = [a - 2 for a in X]
Y = [a - 1 for a in Y]
print X, np.shape(X)
print Y, np.shape(Y)
mlon = np.linspace(-180,180,72)[np.array(X)]
mlat = np.linspace(-90,90,46)[np.array(Y)]
print mlon, mlat

from matplotlib import pyplot as plt
import matplotlib as mpl
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.array([0,0.001,0.01,0.1,0.2,0.4,1.0])
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fignum = 6
if fignum in ptm:
    pl.figure(fignum, figsize=(13,9))
    bm = Basemap(projection = 'cyl', llcrnrlon=min(mlon),llcrnrlat=min(mlat),
                 urcrnrlon=max(mlon),urcrnrlat=max(mlat))
    #lons, lats = np.meshgrid(lon, lat)
    xx, yy = bm(mlon, mlat)
    bm.drawcoastlines(linewidth=1.25, color='k')
    bm.drawparallels(np.array([10,30,50,70]),
                     labels=[1,0,0,0],fontsize=20)
    bm.drawmeridians([-150,-100,-50,0,50,100,150],labels=[0,0,1,0],
                     fontsize=20)
    
    print np.shape(pfnaY)
    bm.pcolor(xx,yy,pfnaY.reshape((len(mlon),len(mlat))).T,cmap=cmap, norm=norm)
    ax = pl.gca()
    divider = mal(ax)
    cax = divider.append_axes('bottom',size='10%', pad = 0.1)
    pl.clim(0.,0.55)
    cbar = pl.colorbar(orientation='horizontal', cax=cax,
                       ticks=[0.,0.001,0.01,.10,.20,.40,1.])
    cbar.ax.set_xticklabels([0.,0.001,0.01,.10,.20,.40,1.],fontsize=20)
    pdump({'lat':mlat,'lon':mlon,'Y':pfnaY.reshape((len(mlon),len(mlat))).T},
          'pfna_data%s.pkl'%tag)
    pl.savefig('pfna_map%s.png'%tag)
    
fignum = 7
if fignum in ptm:
    pl.figure(fignum, figsize=(13,9))
    bm = Basemap(projection = 'cyl', llcrnrlon=min(mlon),llcrnrlat=min(mlat),
                 urcrnrlon=max(mlon),urcrnrlat=max(mlat))
    #lons, lats = np.meshgrid(lon, lat)
    xx, yy = bm(mlon, mlat)
    bm.drawcoastlines(linewidth=1.25, color='k')
    bm.drawparallels(np.array([10,30,50,70]),
                     labels=[1,0,0,0], fontsize=20)
    bm.drawmeridians([-150,-100,-50,0,50,100,150],labels=[0,0,1,0],
                     fontsize=20)    
    print np.shape(pfnaY)
    bm.pcolor(xx,yy,pfoaY.reshape((len(mlon),len(mlat))).T,cmap=cmap,norm=norm)
    ax = pl.gca()
    divider = mal(ax)
    cax = divider.append_axes('bottom',size='10%', pad = 0.1)
    pl.clim(0.,0.55)
    cbar = pl.colorbar(orientation='horizontal', cax=cax,
                       ticks = [0.,0.001,0.01,.10,.20,.40,1.])
    cbar.ax.set_xticklabels([0.,0.001,0.01,.10,.20,.40,1.],fontsize=20)
    pdump({'lat':mlat,'lon':mlon,'Y':pfoaY.reshape((len(mlon),len(mlat))).T},
          'pfoa_data%s.pkl'%tag)
    pl.savefig('pfoa_map%s.png'%tag)

pl.show()
