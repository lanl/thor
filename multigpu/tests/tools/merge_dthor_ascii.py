#!/usr/bin/env python
##
## Merges ASCII TT output from several ranks into a single file
##
import sys, os, glob
import numpy as np

def merge_files(fprefix):
    files = glob.glob(fprefix + "_*.dat")
    nfiles = len(files)
    if nfiles == 0:
        sys.stderr.write("ERROR: couldn't find files\n")
        sys.exit(-2)
    
    if nfiles == 1:
        sys.stderr.write("ERROR: only one file, nothing to merge\n")
        sys.exit(-3)
    
    # check if files are there
    for i in range(nfiles):
        fname = fprefix + "_" + str(i) + ".dat"
        if not os.path.exists(fname):
            sys.stderr.write("ERROR: file "+fname+"  does not exist\n")
            sys.exit(-3)

    # read file sizes and allocate arrays
    fname = fprefix + "_0.dat"
    is_hstack = True
    with open(fname,'r') as f:
        fline = f.readline()[:-1]
        fields = fline.split(" ")
        if (fields[1] == "vstack"):
            is_hstack = False
        elif (fields[1] == "hstack"):
            is_hstack = True
        else:
            sys.stderr.write("ERROR: vstack or hstack are not specified")

        fline = f.readline()[:-1]
        nn = np.array(list(map(int, fline.split(','))),dtype=int)
        fline = f.readline()[:-1]
        rr = np.array(list(map(int, fline.split(','))),dtype=int)
        print (f"tensor modes: {nn}")
        print (f"tensor ranks: {rr}")
    
    d = len(nn)
    cores = []
    for i in range(d):
        cores.append(np.zeros(rr[i]*nn[i]*rr[i+1], dtype=float))
    
    # read the cores
    coffs = np.zeros(d, dtype=int) # cores offsets
    for n in range(nfiles):
        fname = fprefix + "_" + str(n) + ".dat"
        nl = nn
        with open(fname,'r') as f:
            fline = f.readline() # skip first four lines
            fields = fline.split(" ")
            if (fields[1] == "vstack"):
                is_hstack = False
            elif (fields[1] == "hstack"):
                is_hstack = True
            else:
                sys.stderr.write("ERROR: unknown stacking in " + fname + "\n")
                sys.exit(-4)
            fline = f.readline()
            fline = f.readline()
            fline = f.readline() # skip blocking dimension: assume no shredding
            fline = f.readline()[:-1]
            nl = np.array(list(map(int, fline.split(','))),dtype=int)
        vals = np.loadtxt(fname,skiprows=5)
        offset = 0
        for k in range(d):
            clen = rr[k]*nl[k]*rr[k+1] # core length
            if (is_hstack):
                h1 = vals[offset:offset+clen]
            else:
                vvv = np.reshape(vals[offset:offset+clen],[rr[k],nl[k],rr[k+1]],order='F')
                hhh = vvv.transpose((0,2,1))
                #print (vvv.shape, '->', hhh.shape)
                h1  = np.reshape(hhh, [clen], order='F')
            cores[k][coffs[k]:coffs[k]+clen] = h1[0:clen]
            offset += clen
            coffs[k] += clen

    # write output as vstack
    foutname = fprefix + ".dat"
    nlines = 0
    with open(foutname,'w') as f:
        f.write("# vstack\n")

        for i in range(d-1):
            f.write("%d," % nn[i])
            nlines += 1
        f.write("%d\n" % nn[d-1])
        for i in range(d):
            f.write("%d," % rr[i])
            nlines += 1
        f.write("%d\n" % rr[d])
        for k in range(d):
            hhh = np.reshape(cores[k],[rr[k],rr[k+1],nn[k]],order='F')
            vvv = hhh.transpose((0,2,1))
            # print (vvv.shape)
            clen= rr[k]*rr[k+1]*nn[k]
            v1  = np.reshape(vvv, [clen], order='F')
            for i in range(clen):
                f.write("%22.14e\n"% v1[i])
                nlines += 1
    print (f"Output: {foutname}, written lines: {nlines}")


if __name__=='__main__':
    if (len(sys.argv) != 2):
        sys.stderr.write("Merges ASCII dthor output from several ranks into a single file\n");
        sys.stderr.write("Usage: %s <file-prefix>\n"% sys.argv[0])
        sys.exit(-1)

    fnamprefix = sys.argv[1]
    merge_files(fnamprefix)
