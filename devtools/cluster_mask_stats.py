'''
Tool to analyse a .klg file and produce statistics on the size of cluster masks.
'''

from pylab import *

import sys
import re

def compute_stats(fname, itype_filter=None, title_addition="all"):
    logfile = open(fname, 'r').read()
    pattern = r'Cluster mask: cluster (\d+) unmasked (\d+) iterations (\d+)/(\d+) init type (\d+).'
    cluster_mask_lines = re.findall(pattern, logfile)
    cluster, unmasked, localiter, globaliter, itype = zip(*cluster_mask_lines)
    cluster = array(map(int, cluster))
    unmasked = array(map(int, unmasked))
    localiter = array(map(int, localiter))
    globaliter = array(map(int, globaliter))
    itype = array(map(int, itype))
    
    if itype_filter is not None:
        I = itype==itype_filter
        if sum(I)==0:
            return True
        globaliter = globaliter[I]
        # relabel main iterations
        globaliter = cumsum(hstack((0, minimum(diff(globaliter), 1))))
        cluster = cluster[I]
        unmasked = unmasked[I]
        localiter = localiter[I]
        itype = itype[I]
    
    numiter = amax(globaliter)+1
    maxunmasked = amax(unmasked)+1
    allmean = []
    allstd = []
    allmin = []
    allmax = []
    allhist = []
    iters = arange(numiter)
    for iter in iters:
        I = globaliter==iter
        U = unmasked[I]
        allmean.append(mean(U))
        allstd.append(std(U))
        allmin.append(amin(U))
        allmax.append(amax(U))
        b = bincount(U, minlength=maxunmasked)[:maxunmasked]
        b = b*1.0/amax(b)
        allhist.append(b)
    allmean = array(allmean)
    allstd = array(allstd)
    allmin = array(allmin)
    allmax = array(allmax)
    allcummin = minimum.accumulate(allmin)
    allcummax = maximum.accumulate(allmax)
    allhist = array(allhist)
    print allhist.shape
    
    figure(figsize=(10, 8))

    def barhist(U):
        n = amax(unmasked)+1
        b = bincount(U, minlength=n)[:n]
        bar(arange(n)-0.5, b, width=1.0, color='k')
    
    subplot(311)
    #hist(unmasked)
    barhist(unmasked)
    xlabel("Number of unmasked features")
    yticks([])
    title("Frequency over entire run (%s)" % title_addition)
    axis('tight')
    
    subplot(334)
    barhist(unmasked[globaliter==0])
    title('Start')
    yticks([])
    
    subplot(335)
    barhist(unmasked[globaliter==numiter/2])
    title('Middle')
    yticks([])
    
    subplot(336)
    barhist(unmasked[globaliter==numiter-1])
    title('End')
    yticks([])
    
    subplot(313)
    imshow(-allhist.T, interpolation='nearest', aspect='auto', origin='lower left')
    gray()
    fill_between(iters, allmin, allmax, color='b', alpha=0.2)
    plot(iters, allmean, '-b', label='Mean')
    plot(iters, allcummin, '-g', label='Cumulative min')
    plot(iters, allcummax, '-r', label='Cumulative max')
    title("Number of unmasked features over time")
    xlabel("Global iteration number")
    legend(loc='best', fontsize='x-small')
    axis('tight')
    
    tight_layout()
    
    return False
    

if __name__=='__main__':
    
    # DEBUG mode, TODO: remove
    #sys.argv.append('../temp/testsmallish.klg.4')
    #sys.argv.append('../temp/20141008_10000.klg.1')
    
    if len(sys.argv)<=1:
        print 'Usage: cluster_mask_stats filename.klg'
        exit()
    
    fname = sys.argv[1]
    
    compute_stats(fname, None, "all iteration types")
    compute_stats(fname, 0, "main iterations only")
    compute_stats(fname, 1, "K3 iterations only")
    compute_stats(fname, 2, "K2 iterations only")

    show()
    