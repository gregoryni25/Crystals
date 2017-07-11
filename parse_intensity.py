import numpy as np
import matplotlib.pyplot as pl
import time

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

def theoretical_bins(edges, z=None, n=1.82, c=299792458, nn=100, NN=1):
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

def parse_segs_data(D, aline):
    """
    a utility used by parse_segs. gathers relevant line data
    """
    aline = aline.split()  #,'%s','delimiter','\t');

    # hit face 3, the front detector
    if int(aline[4]) == 3:
        return ["front", int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), float(aline[12])]

    # hit face 4, the back detector
    if int(aline[4]) == 4:
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

def parse_segs(fid, numseg):
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

        D.append(parse_segs_data(D, aline))

        if D[idx][0] == "front":
            R.append((1, D[idx][3], pathlength(D, D[idx][2], D[idx][1])))
        if D[idx][0] == "back":
            R.append((0, D[idx][3], pathlength(D, D[idx][2], D[idx][1])))
    
    return R

def parse(filename, numrays):
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
        A += parse_segs(fid, numseg)
 
        idx += 1
    fid.close()

    return A

# returns (number of rays, number hit front, edges, tops)
def specs(parsedarr, plhistarr):
    return (len(parsedarr), sum(parsedarr[:,0]), plhistarr[1], plhistarr[0])
    
"""
total rays : len(parsedarr)
total intensity : sum(parsedarr[:,1])
number front : sum(parsedarr[:,0])
number back : len(parsedarr) - sum(parsedarr[:,0])
edge list : plhistarr[1]
top list : plhistarr[0]

"""    
    
t = time.time()

A = np.array(parse('5000-4.txt',5000))
#B = np.array(parse('5000-3.txt',5000))
#C = np.array(parse('5000-2.txt',5000))
_c_ = 299792458*1000
_n_ = 1.82
_mc_ = _n_/_c_
A2 = pl.hist((A[:,2]*_mc_), bins=50, weights=A[:,1], alpha=0.5, label="four")
#B2 = pl.hist(B[:,2]*_mc_, bins=50, weights=B[:,1], alpha=0.5, label="three")
#C2 = pl.hist(C[:,2]*_mc_, bins=50, weights=C[:,1], alpha=0.5, label="two")

pl.legend()
pl.xlabel("Time (seconds)")
pl.ylabel("Total Intensity")

A3 = specs(A,A2)
#B3 = specs(B,B2)
#C3 = specs(C,C2)


print(time.time()-t)
