import numpy as np

lon4x5=[-182.500, -177.500, -172.500, -167.500, -162.500, -157.500,
        -152.500, -147.500, -142.500, -137.500, -132.500, -127.500,
        -122.500, -117.500, -112.500, -107.500, -102.500,  -97.500,
        -92.500,  -87.500,  -82.500,  -77.500,  -72.500,  -67.500,
        -62.500,  -57.500,  -52.500,  -47.500,  -42.500,  -37.500,
        -32.500,  -27.500,  -22.500,  -17.500,  -12.500,   -7.500,
        -2.500,    2.500,    7.500 ,  12.500,   17.500 ,  22.500,
        27.500,   32.500,   37.500,   42.500,   47.500,   52.500,
        57.500,   62.500,   67.500,   72.500,   77.500,   82.500,
        87.500,   92.500,   97.500,  102.500,  107.500,  112.500,
        117.500,  122.500,  127.500,  132.500, 137.500,  142.500 ,
        147.500,  152.500,  157.500,  162.500,  167.500,  172.500,
        177.500]

lat4x5 = [-90.000, - 88.000, - 84.000, - 80.000, - 76.000, - 72.000,
          - 68.000, - 64.000, -60.000, - 56.000, - 52.000, - 48.000,
          - 44.000, - 40.000, - 36.000, - 32.000, -28.000, - 24.000,
          - 20.000, - 16.000, - 12.000, - 8.000, - 4.000, 0.000,
          4.000, 8.000, 12.000, 16.000, 20.000, 24.000, 28.000,
          32.000, 36.000, 40.000, 44.000, 48.000, 52.000, 56.000,
          60.000, 64.000, 68.000, 72.000, 76.000, 80.000, 84.000,
          88.000, 90.000]

def getGCbox(lat,lon,grid='4x5'):
    """Give i,j GC grid location for given lat, lon."""
    if grid=='4x5':
        for gi,glat in enumerate(lat4x5):
            if glat < lat <= lat4x5[gi+1]:
                i = gi
                break
            else:
                i = len(lat4x5)
        for gj,glon in enumerate(lon4x5):
            if glon < lon <= lon4x5[gj+1]:
                j = gj
                break
            else:
                j = len(lon4x5)
    return i,j

def parse_hysplit(filename):
    """extract lats, lons, heights, hours"""
    data = np.genfromtxt(filename)
    return data[:,9],data[:,10],data[:,11],data[:,8]

def getGCbox_fromfile(filename):
    lats,lons,heights,hours = parse_hysplit(filename)
    Is,Js = [],[]
    for lat,lon in zip(lats,lons):
        #print(lat,lon)
        i,j = getGCbox(lat,lon)
        #print (i,j)
        Is.append(i)
        Js.append(j)
    return Is,Js

if __name__=='__main__':
    ii,jj = getGCbox_fromfile('hysplit_traj.txt')
    print(min(ii),max(ii),min(jj),max(jj))