from pytools import pload
import pylab as pl
import os
import sys
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler



rundir = sys.argv[1]
ensemblename = 'ensemble_summer_highres'

stack = pload(os.path.join(rundir, '%s_data_stack.pkl'%ensemblename))
names = pload(os.path.join(rundir, '%s_data_names.pkl'%ensemblename))
print names

for i in range(13):
    print names[i], i, min(stack[i,:]),max(stack[i,:])


lon = stack[0,:]
lat = stack[1,:]
plotx = stack[3,:]
ploty = stack[2,:]
print lon,lat
x3 = stack[8,:]
x4 = stack[9,:]
x5 = stack[10,:]

# TRANSPOSED FOR DBSCAN!!!
ministack = np.vstack((x3,x4,x5)).T
stdstack = StandardScaler().fit_transform(ministack)
clustered = DBSCAN(eps=0.3, min_samples=10).fit(stdstack)#0.4

core_samples_mask = np.zeros_like(clustered.labels_, dtype=bool)
core_samples_mask[clustered.core_sample_indices_] = True
labels = clustered.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters: %d' % n_clusters_)

import matplotlib.pyplot as plt

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = plt.cm.gist_rainbow(np.linspace(0, 1, len(unique_labels)))
colors = ['b','r','k']

plt.figure(figsize=(12,9))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    XY = np.vstack((plotx,ploty)).T
    xy = XY[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14, alpha=0.6)

    xy = XY[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6, alpha=0.6)
    plt.xlabel('Formation time (days)', fontsize=25)
    plt.ylabel('PFNA Yield', fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

#plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.savefig('clusterspace.png')

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(2)
ax = fig.gca(projection='3d')
#plt.figure(2)
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'
    class_member_mask = (labels == k)

    XYZ = stdstack
    xyz = XYZ[class_member_mask & core_samples_mask]
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:,2],'o', c=col,s=14)

    xyz = XYZ[class_member_mask & ~core_samples_mask]
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:,2],'o', c=col,s=6)
    ax.set_xlabel('log[OH]')
    ax.set_ylabel('log[HO2]')
    ax.set_zlabel('log[NO]')


plt.figure(3,figsize=(18,12))
from mpl_toolkits.basemap import Basemap

mlon = np.linspace(-180,175, 72)[np.array([int(L) for L in lon])]
mlat = np.linspace(-90,90,46)
mlat[-1]=89.
mlat = mlat[np.array([int(L) for L in lat])]

XY = np.vstack((mlon,mlat)).T

bm = Basemap(projection = 'cyl', llcrnrlon=min(mlon),llcrnrlat=min(mlat),
             urcrnrlon=max(mlon),urcrnrlat=max(mlat))
#lons, lats = np.meshgrid(lon, lat)
x, y = bm(mlon, mlat)

colorarray = np.zeros((len(x),len(y)))

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'
        siz = 6
    else:
        siz = 11
    class_member_mask = (labels == k)

    XY = np.array([x,y]).T
    xy = XY[class_member_mask & core_samples_mask]
    bm.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
            markersize=siz, alpha=0.7)

    xy = XY[class_member_mask & ~core_samples_mask]
    bm.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
            markersize=siz, alpha=0.7)

bm.drawcoastlines(linewidth=2, color='k')
bm.drawparallels(np.array([10,30,50,70]),
                 labels=[1,0,0,0], fontsize=20)
bm.drawmeridians(np.arange(min(mlon),max(mlon),60),labels=[0,0,0,1],fontsize=20)
plt.savefig('mapclusters.png')

plt.show()
