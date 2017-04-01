# -*- coding: utf-8 -*-
"""
Created on Sat May 11 18:06:23 2013

@author: thackray
"""

import os

from csv import DictReader
import numpy as np
import _pickle as pickle

try:
    import pylab as pl
    PL = True
except ImportError:
    PL = False
try:
    from mpl_toolkits.basemap import Basemap
    BASEMAP = True
except ImportError:
    BASEMAP = False
try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable as mal
    MAL = True
except ImportError:
    MAL = False

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',
          'Sep', 'Oct', 'Nov', 'Dec']

def compare(x,y, eps=1e-6):
    return abs(x-y) < eps

def savefig(figname):
    """Takes a .png-less fig name and saves it with the .png extension and 
    redued whitespace around the meat of the figure"""
    pl.savefig(figname+'.png', bbox_inches = 'tight')


def pload(filename, asarray = False,encoding="ASCII"):
    """Wrapper to make un-pickling friendlier.

    Arguments:
    filename: string file name to open

    Keyword Arguments:
    asarray: (Default False) will return pickled object as an array if possible

    Returns:
    Object that was pickled
    """

    f = open(filename, 'rb')
    d = pickle.load(f,encoding=encoding)
    
    if asarray:
        return np.array(d)
    else:
        return d

def pdump(obj, filename):
    """Wrapper to make pickling friendlier.

    Arguments:
    obj: object to be pickled
    filename: name of file to pickle object under

    Returns:
    Nothing
    """
    f = open(filename, 'wb')
    pickle.dump(obj,f)
    return

def listget(listy, i, default=0):
    try:
        return listy[i]
    except IndexError:
        return default

def dict_to_csv(dictionary, csvname):
    """Save a dictionary to csv format.

    Arguments:
    dictionary: dictionary to unpack into csv
    csvname: string filename

    Returns:
    None
    """

    fieldnames = dictionary.keys()
    fieldlengths = [len(dictionary[field]) for field in fieldnames]
    with open(csvname, 'w') as f:
        f.write(','.join(fieldnames)+'\n')
        for i in range(max(fieldlengths)):
            line = [str(listget(dictionary[fld],i,'')) for fld in fieldnames]
            f.write(','.join(line)+'\n')
    return


def dictcsv(csvname, fieldnames = None, arrays = False):
    """Reading csv files into a dictionary.

    Arguments:
    csvname: string filename
    
    Keyword Arguments:
    fieldnames: list of csv column names.  If none, first column of the file
                being read will be used.
    arrays: Whether or not to return csv contents as a dict of arrays

    Returns: 
    dictionary of columns as numpy arrays, keys are fieldnames    
    """

    fileobj = open(csvname, 'rU')
    DR = DictReader(fileobj, fieldnames = fieldnames)
    
    fields = DR.fieldnames
    l = DR.next()
    dicty = {}
    for f in fields:
        try:
            dicty[f] = [float(l[f])]
        except (TypeError, ValueError):
            dicty[f] = [l[f]]
    for row in DR:
        for f in fields:
            try:
                dicty[f].append(float(row[f]))
            except (TypeError, ValueError):
                dicty[f].append(row[f])
    if arrays:
        for key in dicty:
            dicty[key] = np.array(dicty[key])
            
    return dicty

if PL and BASEMAP:
    ### GEOSMAP in GEOSPLOT does this better

    def plot_map(lat, lon, data, normalized=False):
        """Plots data on a Basemap.
        
        Arguments:
        lat: array of lat locations of data
        lon: array of lon locations of other axis of data
        data: 2D lat-lon grid of data to put on map
        
        Keyword Arguments:
        normalized: if True, show a 0:1 map of ratios of max value
        
        Returns:
        Nothing
        """

        bm = Basemap(projection = 'mill', llcrnrlon=min(lon),
                     llcrnrlat=min(lat), urcrnrlon=max(lon),
                     urcrnrlat=max(lat))
        lons, lats = np.meshgrid(lon, lat)
        x, y = bm(lons, lats)

        pl.figure(figsize=(18,14))
        ax = pl.gca()
        bm.drawcoastlines(linewidth=1.25, color='white')
        bm.drawparallels(np.array([-70,-50,-30,-10,10,30,50,70]),
                         labels=[1,0,0,0])
        bm.drawmeridians(np.arange(min(lon),max(lon),60),labels=[0,0,0,1])
        
        if normalized:
            bm.pcolor(x,y,data.T/np.abs(np.max(data)))
        else:
            bm.pcolor(x,y,data.T)

        if MAL:
            divider = mal(ax)
            cax = divider.append_axes("right", size='5%', pad=0.05)
            pl.colorbar(cax=cax)
        else:
            pl.colorbar()

        return

if __name__=='__main__':
    dict_to_csv({'a':range(6),'b':range(3)},'test_write.csv')
