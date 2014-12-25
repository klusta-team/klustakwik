'''
Usage:

python data_mask_stats.py name.fmask.n features_per_channel max_difference

Prints the number of unique binary masks in the file.
'''
from numpy import *
import sys
import time
from txt_tools import *

def unique_rows(masks):
    # masks.shape = (numpoints, numfeatures)
    I = lexsort(masks.T)
    masks = masks[I]
    d = diff(masks, axis=0)
    d2, = (sum(abs(d), axis=1)>0).nonzero()
    d3 = hstack((0, 1+d2))
    masks = masks[d3]
    return masks

def num_unique_rows(arr):
    return unique_rows(arr).shape[0]

def greedy_reduce(masks, candidates, max_difference):
    # stage 1: extend candidates set
    print 'Searching for new candidates'
    n_masks = sum(masks, axis=1)
    n_candidates = sum(candidates, axis=1)    
    new_candidates = [candidates]
    for i in xrange(masks.shape[0]):
        print 'Searching for new candidates', i,
        nfound = 0
        A = masks[i, :]
        ND = n_masks[i]-n_candidates
        J, = (abs(ND)<=1).nonzero()
        MI = A[newaxis, :]
        CJ = candidates[J, :]
        C = MI | CJ
        NC, = ((sum(C & -MI, axis=1)<=1) & (sum(C & -CJ, axis=1)<=1)).nonzero()
        new_candidates.append(C[NC, :])
        print 'found', len(NC)
    new_candidates = vstack(new_candidates)    
    candidates = unique_rows(candidates)
    # stage 2: compute coverings
    n_candidates = sum(candidates, axis=1)
    covering = [set([]) for _ in xrange(n_candidates.shape[0])]
    covered_by = [set([]) for _ in xrange(n_masks.shape[0])]
    for i in xrange(masks.shape[0]):
        print 'Computing coverings for', i,
        A = masks[i, :]
        nfound = 0
        ND = n_candidates-n_masks[i]
        J, = (ND<=max_difference).nonzero()
        MI = A[newaxis, :]
        CJ = candidates[J, :]
        C = MI <= CJ
        NC, = C.all(axis=1).nonzero()
        for j in J[NC]:
            covering[j].add(i)
            covered_by[i].add(j)
            nfound += 1
        print 'found', nfound
    # stage 3: greedily select best coverings
    remaining = set(arange(masks.shape[0]))
    reduced = []
    num_covered = [len(C) for C in covering]
    while remaining:
        print 'Assigning coverings, number of remaining masks', len(remaining)
        k = argmax(num_covered)
        remaining.difference_update(covering[k])
        reduced.append(candidates[k, :])
        for i in covering[k].copy():
            for j in covered_by[i]:
                assert i in covering[j]
                covering[j].remove(i)
                num_covered[j] -= 1
    return array(reduced)
    
def num_reduced_masks(masks, max_difference):
    masks = unique_rows(masks)
    candidates = masks
    all_nrm = []
    for k in xrange(1, max_difference+1):
        candidates = greedy_reduce(masks, candidates, k)
        all_nrm.append(candidates.shape[0]*1)
    return all_nrm

def analyse_fmask(fname, FPC, max_difference):
    fmask = loadtxt(fname, skiprows=1, dtype=float)
    fmask = fmask[:, ::FPC]
    numspikes, numfeatures = fmask.shape
    masks = fmask>0
    n = num_unique_rows(masks)
    nr = num_reduced_masks(masks, max_difference) 
    print 'Number of unique binary masks is %d (%.2f %% of possible)' % (n, n*100.0/numspikes)
    print 'Number of channels is %d, number of points is %d' % (fmask.shape[1], fmask.shape[0])
    for k, nrk in enumerate(nr):
        print 'Number of unique reduced binary masks with max difference %d is %d (%.2f %% of possible)' % (k+1, nrk, nrk*100.0/numspikes)

if __name__=='__main__':
    
    # DEBUG mode, TODO: remove
    #sys.argv.extend(['../temp/testsmallish.fmask.4', '3', '1'])
    #sys.argv.extend(['../temp/20141008_10000.fmask.1', '3', '4'])

    if len(sys.argv)==4:
        fname = sys.argv[1]
        FPC = int(sys.argv[2])
        max_difference = int(sys.argv[3])
        analyse_fmask(fname, FPC, max_difference)
        exit()
    
    print __doc__
