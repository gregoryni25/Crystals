import numpy as np
import matplotlib.pyplot as pl
import time

def edge2cen(edges):
    return (edges[:-1] + edges[1:])/2

def dens_m(z, n, c, t):
    '''density fcn for unit power'''
    coeff = z*n/c
    out = coeff/2 / t**2
    out[t < coeff] = 0
    return out

def trape_dens(edge1, edge2, z, n, c, nn):
    t = np.linspace(edge1, edge2, nn)
    y = dens_m(z, n, c, t)
    return np.trapz(y, t)

def theoretical_bins(edges, z=None, n=1.82, c=299792458000, nn=100, NN=10000):
    '''
    edges       bin edges in time (s)
    z           distance from source to end cap
    n           refractive index
    c           speed of light
    nn          number of points used in trapezoidal integration
    NN          total number of rays
    output:
    contents    bin contents
    '''
    if z is None:
        z = edges[0]*c/n
    contents = np.empty(len(edges)-1)
    for i in range(len(contents)):
        contents[i] = trape_dens(edges[i], edges[i+1], z, n, c, nn)
    return NN * contents

def pathlength(D, fcoor, parent):
    """
    takes in temp Data array, first coordinate and parent, and returns the
    pathlength from root ray (ray 0) to that ray
    """
    coor = D[parent][1]
    
    # base case: at segment 0
    if D[parent][0] == "base":
        return np.linalg.norm(np.subtract(fcoor, coor))
    
    # recursive case: if D[parent][0] is positive, than it points to
    # its parent, otherwise it contains pathlength already
    if D[parent][0] >= 0:
        D[parent][0] = -pathlength(D, coor, D[parent][0])
    
    return np.linalg.norm(np.subtract(fcoor, coor)) - D[parent][0]    

def parse_segs_data(D, aline, level=False):
    """
    a utility used by parse_segs. gathers relevant line data
    """
    aline = aline.split()  #,'%s','delimiter','\t');

    # hit face 3, the front detector
    if int(aline[4]) == 3:
        if level:
            return ["front", int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), int(aline[2])]
        return ["front", int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), float(aline[12])]

    # hit face 4, the back detector
    if int(aline[4]) == 4:
        if level:
            return ["back", int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), int(aline[2])]
        return ["back", int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), float(aline[12])]

    # if it didn't hit detector 3 or 4 but still terminated, then don't store
    if aline[6] == '*---':
        return [0]
    
    # base case
    if aline[0] == '0':
        return ["base", (float(aline[9]), float(aline[10]),
             float(aline[11])), 0]

    return [int(aline[1]), (float(aline[9]), float(aline[10]),
             float(aline[11]))]

def parse_segs(fid, numseg, level=False):
    """
    a utility used by parse. goes through a ray and finds relevant data
    for rays that hit detector 3 or 4
    """
    # temporary storage array
    D = []

    # answer array
    R = []

    for idx in range(numseg+1):
        aline = fid.readline()

        D.append(parse_segs_data(D, aline, level))

        if D[idx][0] == "front":
            R.append((1, D[idx][3], pathlength(D, D[idx][2], D[idx][1])))
        if D[idx][0] == "back":
            R.append((0, D[idx][3], pathlength(D, D[idx][2], D[idx][1])))
    
    return R

def parse(filename, numrays, level=False):
    """
    Reads file, returns an array of rays' (front/back, level, pathlength)
    REQUIRES: face 4 is back detector, face 3 is front detector
    returns 0 if hit face 4, 1 if hit face 3
    """
    # filename='crystal.1.txt'
    #numrays=1000;
    fid = open(filename, 'r', encoding='utf-16-le')  # 539294 19328874 643694
    # remove 12 lines
    for idx in range(12):
        fid.readline()

    # initialize answer array
    A = []

    idx = 0
    while idx < numrays:  # for idx in range(numrays):
        # should be able to read for numrays times

        # find numseg
        aline = fid.readline()  # 3411771

        if aline == '\n':  # don't process, don't count, if line is empty
            continue
        elif aline == '':  # eof 10140281
            if idx < numrays:
                print('in parse. eof while numrays-idx=', numrays-idx)
            break

        numseg = int(aline.split()[6])
        
        # head of table for each ray
        fid.readline()
        
        # go through segs, accumulate pathlength
        # this results in [[r1],[r2],[r3]...]
        # A.append(parse_segs(fid, numseg))
        # this results in [r1,r2,r3...]
        A += parse_segs(fid, numseg, level)
 
        idx += 1
    fid.close()

    return A

# returns (number of rays, number hit front, edges, tops)
def specs(parsedarr, plhistarr):
    return (len(parsedarr), sum(parsedarr[:,0]), plhistarr[1], plhistarr[0])

#returns an array with only rays that hit the front detector
def onlyfront(parsedarr):
    indexarray = []
    for idx in range(len(parsedarr)):
        if parsedarr[idx][0] == 0:
            indexarray.append(idx)
    return np.delete(parsedarr, indexarray, 0)

#returns an array with only rays that hit the back detector
def onlyback(parsedarr):
    indexarray = []
    for idx in range(len(parsedarr)):
        if parsedarr[idx][0] == 1:
            indexarray.append(idx)
    return np.delete(parsedarr, indexarray, 0)
            

"""
total rays : len(parsedarr)
total intensity : sum(parsedarr[:,1])
number front : sum(parsedarr[:,0])
number back : len(parsedarr) - sum(parsedarr[:,0])
edge list : plhistarr[1]
top list : plhistarr[0]

"""    
    
t = time.time()

data_path = '../Crystal Data/'

_numrays_ = 5000
A = np.array(parse('%d-2-2.txt' % _numrays_ ,_numrays_, level=True))
#B = np.array(parse(data_path+'5000-3.txt',5000))
#C = np.array(parse(data_path+'5000-2.txt',5000))
_c_ = 299792458*1000 #since distance measured in millimeters
_n_ = 1.82
_mc_ = _n_/_c_

iflog = True

Af = onlyfront(A)
Ab = onlyback(A)

_split_=2
_levelweight_=(1/_split_)
#1/split ^ level
#, weights=_levelweight_ ** (Af[:,1]-1)

#, weights=Af[:,1]
#, range=(np.amin(Af[:,2]*_mc_), 2.5e-9)
A1 = pl.hist((Af[:,2]*_mc_), bins=50, range=(np.amin(Af[:,2]*_mc_), 2.5e-9), weights=_levelweight_ ** (Af[:,1]-1), alpha=0.5, label="front",
             log=iflog)

A2 = pl.hist((Ab[:,2]*_mc_), bins=50, range=(np.amin(Ab[:,2]*_mc_), 2.5e-9), weights=_levelweight_ ** (Ab[:,1]-1), alpha=0.5, label="back",
             log=iflog)

#B2 = pl.hist(B[:,2]*_mc_, bins=50, weights=B[:,1], alpha=0.5, label="3",
#             log=iflog)
#C2 = pl.hist(C[:,2]*_mc_, bins=50, weights=C[:,1], alpha=0.5, label="2",
#             log=iflog)


pl.legend()
pl.xlabel("Time (seconds)")
pl.ylabel("Total Intensity")

'''
bin_edges = A1[1]
A_th1 = theoretical_bins(bin_edges,NN=_numrays_)
pl.plot(edge2cen(bin_edges), A_th1, '.')
bin_edges = A2[1]
A_th2 = theoretical_bins(bin_edges,NN=_numrays_)
pl.plot(edge2cen(bin_edges), A_th2, '.')
'''
#sumAth = A_th1 + A_th2

#pl.plot(edge2cen(bin_edges), sumAth, '.')
if iflog:
    pl.gca().set_yscale('log')
    
print(time.time()-t)
