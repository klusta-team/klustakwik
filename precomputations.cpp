/*
 * precomputations.cpp
 *
 * Used for certain precomputations designed to speed up the main loops.
 *
 *  Created on: 19 Nov 2011
 *      Author: dan
 */

#include "klustakwik.h"
#include<algorithm>

using namespace std;

// Handles doing all the Initial precomputations once the data has been loaded into
// the class.
void KK::DoInitialPrecomputations()
{
	if(UseClusterPenalty)
	{
		// Precompute the indices of the unmasked dimensions for each point
		ComputeUnmasked();
		// Compute the order of points to consider that minimises the number of
		// times the mask changes
		ComputeSortIndices();
		// Now compute the points at which the mask changes in sorted order
		ComputeSortedUnmaskedChangePoints();
        // Compute the sum of the masks/float masks for each point (used for computing the cluster penalty)
        PointMaskDimension(); 
	}
	if(UseDistributional)
	{
		// Precompute the noise means and variances
		ComputeNoiseMeansAndVariances();
	    ComputeCorrectionTermsAndReplaceData();
	}
}

// Handles doing all the precomputations once the data has been loaded into
// the class. Note that this function has to be called after Data or Masks
// has changed. This is called by TrySplits() and Cluster()

void KK::DoPrecomputations()
{
	if(UseClusterPenalty)
	{
		// Precompute the indices of the unmasked dimensions for each point
		ComputeUnmasked();
		// Compute the order of points to consider that minimises the number of
		// times the mask changes
		ComputeSortIndices();
		// Now compute the points at which the mask changes in sorted order
		ComputeSortedUnmaskedChangePoints();
	}

	if(UseDistributional)
	{
		ComputeCorrectionTermsAndReplaceData();
	}
}

void KK::ComputeUnmasked()
{
	int i=0;
	if(Unmasked.size() || UnmaskedInd.size())
	{
		Error("Precomputations have already been done, this indicates a bug.\n");
		Error("Error occurred in ComputeUnmasked().\n");
		abort();
	}
	for(int p=0; p<nPoints; p++)
	{
		UnmaskedInd.push_back(i);
		for(int j=0; j<nDims; j++)
		{
			if(Masks[p*nDims+j])
			{
				Unmasked.push_back(j);
				i++;
			}
		}
	}
	UnmaskedInd.push_back(i);
}

// This function computes the points at which the mask changes if we iterate
// through the points in sorted order defined by ComputeSortIndices(). After
// the function is called, SortedMaskChange[SortedIndices[q]] is true if
// the mask for SortedIndices[q] is different from the mask for
// SortedIndices[q-1]
void KK::ComputeSortedUnmaskedChangePoints()
{
	if(SortedMaskChange.size()>0)
	{
		Error("Precomputations have already been done, this indicates a bug.\n");
		Error("Error occurred in ComputeSortedUnmaskedChangePoints().\n");
		abort();
	}
	SortedMaskChange.resize(nPoints);
	SafeArray<int> safeSortedMaskChange(SortedMaskChange, "CSUCP:SMC");
	// The first point when we iterate through the points in sorted order is
	// SortedIndices[0] and we consider the mask as having 'changed' for this
	// first point, because we use the mask having changed to signal that
	// we should recompute the matrices that depend on the masks.
	safeSortedMaskChange[SortedIndices[0]] = true;
	SafeArray<int> oldmask(Masks, SortedIndices[0]*nDims,
			"ComputeSortedUnmaskedChangePoints:oldmask");
	int numchanged = 0;
	for(int q=1; q<nPoints; q++)
	{
		int p = SortedIndices[q];
		SafeArray<int> newmask(Masks, p*nDims,
				"ComputeSortedUnmaskedChangePoints:newmask");
		bool changed = false;
		for(int i=0; i<nDims; i++)
		{
			if(newmask[i]!=oldmask[i])
			{
				oldmask = newmask;
				changed = true;
				numchanged++;
				break;
			}
		}
		safeSortedMaskChange[p] = changed;
	}
}

///////////////// SORTING /////////////////////////////////////////////////
// Comparison class, the operator()(i, j) function is used to provide the
// comparison i<j passed to stl::sort. Here i and j are the indices of two
// points, and i<j means that mask_i < mask_j in lexicographical order (we
// could change the ordering, as long as different masks are not considered
// equal).
class KKSort
{
public:
	KK *kk;
	KKSort(KK *kk) : kk(kk) {};
	bool operator()(const int i, const int j) const;
};

// Less than operator for KK.Masks, it's just a lexicographical comparison
bool KKSort::operator()(const int i, const int j) const
{
	int nDims = kk->nDims;
	for(int k=0; k<nDims; k++)
	{
		int x = kk->Masks[i*nDims+k];
		int y = kk->Masks[j*nDims+k];
		if(x<y) return true;
		if(x>y) return false;
	}
	return false;
}

/*
 * This function computes the order in which indices should be considered in
 * order to minimise the number of times the mask changes. We do this simply
 * by creating an array SortedIndices=[0,1,2,...,nPoints-1] and then sorting
 * this array where i<j if mask_i<mask_j in lexicographical order. The sorting
 * is performed by stl::sort and the comparison function is provided by the
 * KKSort class above.
 *
 * The optional force flag forces a recomputation of the sorted indices, which
 * is necessary only in TrySplits(), but we should probably change this by
 * refactoring.
 */
void KK::ComputeSortIndices()
{
	KKSort kksorter(this);
	if(SortedIndices.size())
	{
		Error("Precomputations have already been done, this indicates a bug.\n");
		Error("Error occurred in ComputeSortIndices().\n");
		abort();
	}
	SortedIndices.resize(nPoints);
	for(int i=0; i<nPoints; i++)
		SortedIndices[i] = i;
	stable_sort(SortedIndices.begin(), SortedIndices.end(), kksorter);
}

//void KK::ComputeNoiseMeansAndVariances()
//{
	//For TrySplits
	//maintain noise mean and variance of each channel
	// maintain number of masked points in each channel
//	Output("ComputeNoiseMeansandVariances ");

//	NoiseMean.resize(nDims);
//	NoiseVariance.resize(nDims);
//	nMasked.resize(nDims);

//}


void KK::ComputeNoiseMeansAndVariances()
{
	// compute noise mean and variance of each channel
	// compute number of masked points in each channel
	Output("Masked EM: Computing Noise Means and Variances \n -----------------------------------------");

	NoiseMean.resize(nDims);
	NoiseVariance.resize(nDims);
	nMasked.resize(nDims);
	
//	for(int i=0; i<nDims; i++)
//	{	NoiseMean[i] = 0;
//		NoiseVariance[i] = 0;
//		nMasked[i] = 0;
//	}
		
	
	for(int p=0; p<nPoints; p++)
		for(int i=0; i<nDims; i++)
			if(!Masks[p*nDims+i])
			{
				scalar thisdata = Data[p*nDims+i];
				NoiseMean[i] += thisdata;
			//	NoiseVariance[i] += thisdata*thisdata; // sum of squares
				nMasked[i]++;
			}
	for(int i=0; i<nDims; i++)
	{
		if(nMasked[i]==0)
		{
			NoiseMean[i] = 0.0;
            NoiseVariance[i] = 0;
	//		NoiseVariance[i] = 1.0;
		} else
		{
			NoiseMean[i] /= (scalar)nMasked[i];
		//	NoiseVariance[i] /= (scalar)nMasked[i]; // E[X^2]
		//	NoiseVariance[i] -= NoiseMean[i]*NoiseMean[i]; // -E[X]^2
		}
	}
	
	
	for(int p=0; p<nPoints; p++)
		for(int i=0; i<nDims; i++)
			if(!Masks[p*nDims+i])
			{	scalar thisdata = Data[p*nDims+i];
				NoiseVariance[i] += (thisdata-NoiseMean[i])*(thisdata-NoiseMean[i]); 
				
			}
	
	for(int i=0; i<nDims; i++)
	{   
        if(nMasked[i]==0)
		{    NoiseVariance[i]= 0;
		    
        }else {
            NoiseVariance[i] /= (scalar)nMasked[i];
        }

	}
	
//	for(int i=0; i<nDims; i++)
//	{	Output(" NoiseMean[%d] = %f",i,NoiseMean[i]);
//		Output(" NoiseVariance[%d] = %f",i,NoiseVariance[i]);
//		Output(" nMasked[%d] = %d",i,nMasked[i]);
		
//	}
}

void KK::ComputeCorrectionTermsAndReplaceData()
{
	for(int p=0; p<nPoints; p++)
		for(int i=0; i<nDims; i++)
		{
			scalar x = Data[p*nDims+i];
			scalar w = FloatMasks[p*nDims+i];
			scalar nu = NoiseMean[i];
			scalar sigma2 = NoiseVariance[i];
			scalar y = w*x+(1-w)*nu;
			scalar z = w*x*x+(1-w)*(nu*nu+sigma2);
			CorrectionTerm[p*nDims+i] = z-y*y;
			Data[p*nDims+i] = y;
		}
}

//SNK PointMaskDimension() computes the sum of the masks/float masks for each point 
void KK::PointMaskDimension() 
{
	int i,p;
	
    
	
    for (p=0; p<nPoints; p++) 
    {
        UnMaskDims[p]=0;
        for (i=0;i<nDims;i++)
        {
            UnMaskDims[p] += FloatMasks[p*nDims+i];
        }
        if (Debug)
        {
            Output("UnMaskDims[%d] = %f ",p,UnMaskDims[p]);
        }
        
    }	
	
}
