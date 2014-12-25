from numpy import *
import sys
import time

__all__ = ['loadtxt', 'savetxt']

try:
    # much faster than numpy's loadtxt and savetxt
    import pandas
    def loadtxt(fname, skiprows=0, dtype=int):
        return pandas.read_csv(fname, sep=' ', skiprows=skiprows, dtype=dtype, header=None).values
    def savetxt(fname, arr, dtype=int, header=None, comments=None):
        df = pandas.DataFrame(arr)
        with open(fname, 'w') as f:
            if header is not None:
                f.write(header+'\n')
            df.to_csv(f, sep=' ', header=False, index=False)
except ImportError:
    # slightly faster than numpy's loadtxt
    def loadtxt(fname, skiprows=0, dtype=int):
        with open(fname, 'r') as f:
            for _ in xrange(skiprows):
                f.readline()
            count = 0
            for line in f:
                if count==0:
                    numdims = len(line.split(' '))
                count += 1
        out = zeros((count, numdims), dtype=dtype)
        with open(fname, 'r') as f:
            for _ in xrange(skiprows):
                f.readline()
            i = 0
            for line in f:
                x = array(line.split(' '), dtype=dtype)
                out[i, :] = x
                i += 1
        return out
