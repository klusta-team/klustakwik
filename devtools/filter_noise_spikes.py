'''Usage:

To remove noise spikes (defined as those spikes which have more than
max_unmasked unmasked channels):

  python filter_noise_spikes.py filename.fet.n max_unmasked
  
You will need to ensure that there is also a filename.fmask.n file. This will
create three new files:

  filename.filtered.fet.n
  filename.filtered.fmask.n
  filename.filtered.removed.n
  
Now run klustakwik on these (e.g. klustakwik filename.filtered n ...). This will generate a file:

  filename.filtered.clu.n
  
The spike numbers on this will not correspond to your original .fet file, so to convert back, run:

  python filter_noise_spikes.py filename.filtered.clu.n
  
This will generate a new file:

  filename.clu.n
  
This file will correctly correspond to your original .fet file.
'''
from numpy import *
import sys


def filter_spikes(fname, max_unmasked):
    fname = fname.replace('.fet.', '?')
    basefname, shanknum = fname.split('?')
    
    fet_fname_in = basefname+'.fet.'+shanknum
    fmask_fname_in = basefname+'.fmask.'+shanknum
    fet_fname_out = basefname+'.filtered.fet.'+shanknum
    fmask_fname_out = basefname+'.filtered.fmask.'+shanknum
    removed_fname_out = basefname+'.filtered.removed.'+shanknum
    
    fet = loadtxt(fet_fname_in, skiprows=1, dtype=int)
    fmask = loadtxt(fmask_fname_in, skiprows=1, dtype=float)
    
    numspikes, numfeatures = fet.shape
        
    sum_fmask = sum(fmask, axis=1)
    cond = sum_fmask<=max_unmasked
    indices_removed, = (sum_fmask>max_unmasked).nonzero()
    num_removed = len(indices_removed)
    num_kept = numspikes-num_removed
    
    savetxt(removed_fname_out, hstack([numspikes, numfeatures, indices_removed]), '%d')
    savetxt(fet_fname_out, fet[cond, :], '%d', header=str(numfeatures), comments='')
    savetxt(fmask_fname_out, fmask[cond, :], '%d', header=str(numfeatures), comments='')
    
    print 'Number of spikes = %d, number of features = %d' % (numspikes, numfeatures)
    print 'Number of spikes removed = %d (%.2f %%)' % (num_removed, num_removed*100.0/numspikes)
    print 'Number of spikes kept = %d (%.2f %%)' % (num_kept, num_kept*100.0/numspikes)

        
def convert_clu(fname):
    fname = fname.replace('.filtered.clu.', '?')
    basefname, shanknum = fname.split('?')
    
    removed_fname_in = basefname+'.filtered.removed.'+shanknum
    clu_fname_in = basefname+'.filtered.clu.'+shanknum
    clu_fname_out = basefname+'.clu.'+shanknum
    fet_fname_in = basefname+'.fet.'+shanknum
    
    removed_indices = loadtxt(removed_fname_in, dtype=int)
    (numspikes, numfeatures), removed_indices = removed_indices[:2], removed_indices[2:]
    clu = loadtxt(clu_fname_in, skiprows=1, dtype=int)

    cond = ones(numspikes, dtype=bool)
    cond[removed_indices] = False

    clu_out = zeros(numspikes)
    clu_out[cond] = clu
    savetxt(clu_fname_out, clu_out, '%d', header=str(amax(clu)), comments='')
        

if __name__=='__main__':
    # DEBUG mode, TODO: remove
    #sys.argv.extend(['../temp/testsmallish.fet.4', '11'])
    #sys.argv.append('../temp/testsmallish.filtered.clu.4')

    if len(sys.argv)==2:
        fname = sys.argv[1]
        if '.filtered.clu.' in fname:
            convert_clu(fname)
            exit()
    
    if len(sys.argv)==3:
        fname = sys.argv[1]
        max_unmasked = float(sys.argv[2])
        if '.fet.' in fname:
            filter_spikes(fname, max_unmasked)
            exit()

    print __doc__
