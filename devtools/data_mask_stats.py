'''
Usage:

python data_mask_stats.py name.fmask.n

Prints the number of unique binary masks in the file.
'''
from numpy import *
import sys
import time
from txt_tools import *

def num_unique_rows(arr):
    order = lexsort(arr.T)
    arr = arr[order]
    d = diff(arr, axis=0)
    d2 = sum(abs(d), axis=1)>0
    return 1+sum(d2)

def analyse_fmask(fname):
    fmask = loadtxt(fname, skiprows=1, dtype=float)
    numspikes, numfeatures = fmask.shape
    masks_nonzero = fmask>0
    masks_one = fmask==1
    n0 = num_unique_rows(masks_nonzero)
    n1 = num_unique_rows(masks_one)
    print 'Number of unique binary masks with fmask>0 is %d (%.2f %% of possible)' % (n0, n0*100.0/numspikes) 
    print 'Number of unique binary masks with fmask==1 is %d (%.2f %% of possible)' % (n1, n1*100.0/numspikes)

if __name__=='__main__':
    
    # DEBUG mode, TODO: remove
    #sys.argv.extend(['../temp/testsmallish.fmask.4'])

    if len(sys.argv)==2:
        fname = sys.argv[1]
        analyse_fmask(fname)
        exit()
    
    print __doc__
