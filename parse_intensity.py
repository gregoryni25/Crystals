import numpy as np
import matplotlib.pyplot as pl
import time

def pathlength(D, fcoor, parent):
    """
    takes in temp Data array, first coordinate and parent, and returns the
    pathlength from root ray (ray 0) to that ray
    """
    coor = D[parent][1]
    return D[parent][2] + np.linalg.norm(np.subtract(fcoor, coor))


def parse_segs_data(D, aline):
    """
    a utility used by parse_segs. gathers relevant line data
    """
    aline = aline.split()  #,'%s','delimiter','\t');

    # hit face 3, the front detector
    if int(aline[4]) == 3:
        return [-3, int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), float(aline[12])]
    # hit face 4, the back detector
    if int(aline[4]) == 4:
        return [-4, int(aline[1]), (float(aline[9]), float(aline[10]),
                 float(aline[11])), float(aline[12])]

    # if it didn't hit detector 3 or 4 but still terminated, then don't store
    if aline[6] == '*---':
        return [0]
    
    # base case
    if aline[0] == '0':
        return [int(aline[1]), (float(aline[9]), float(aline[10]),
             float(aline[11])), 0]
        

    return [int(aline[1]), (float(aline[9]), float(aline[10]),
             float(aline[11])), pathlength(D, (float(aline[9]), 
                   float(aline[10]), float(aline[11])), int(aline[1]))]

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
        if D[idx][0] == -3:
            R.append((1, D[idx][3], pathlength(D, D[idx][2], D[idx][1])))
        if D[idx][0] == -4:
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
t = time.time()

A = np.array(parse('5000-4.txt',5000))
print(time.time()-t)
t = time.time()
B = np.array(parse('5000-3.txt',5000))
print(time.time()-t)
t = time.time()
C = np.array(parse('5000-2.txt',5000))
print(time.time()-t)
pl.hist(A[:,2], bins=20, weights=A[:,1], alpha=0.5, label="four")
pl.hist(B[:,2], bins=20, weights=B[:,1], alpha=0.5, label="three")
pl.hist(C[:,2], bins=20, weights=C[:,1], alpha=0.5, label="two")
pl.legend()
print(time.time()-t)