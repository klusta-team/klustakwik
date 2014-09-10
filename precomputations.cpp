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
    if(UseDistributional)
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
    if(UseDistributional)
    {
        // Precompute the indices of the unmasked dimensions for each point
        ComputeUnmasked();
        // Compute the order of points to consider that minimises the number of
        // times the mask changes
        ComputeSortIndices();
        // Now compute the points at which the mask changes in sorted order
        ComputeSortedUnmaskedChangePoints();
        ComputeCorrectionTermsAndReplaceData();
    }
}

void KK::ComputeUnmasked()
{
    integer i=0;
    if(Unmasked.size() || UnmaskedInd.size())
    {
        Error("Precomputations have already been done, this indicates a bug.\n");
        Error("Error occurred in ComputeUnmasked().\n");
        abort();
    }
    for(integer p=0; p<nPoints; p++)
    {
        UnmaskedInd.push_back(i);
        for(integer j=0; j<nDims; j++)
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
    SafeArray<integer> safeSortedMaskChange(SortedMaskChange, "CSUCP:SMC");
    // The first point when we iterate through the points in sorted order is
    // SortedIndices[0] and we consider the mask as having 'changed' for this
    // first point, because we use the mask having changed to signal that
    // we should recompute the matrices that depend on the masks.
    safeSortedMaskChange[SortedIndices[0]] = true;
    SafeArray<integer> oldmask(Masks, SortedIndices[0]*nDims,
            "ComputeSortedUnmaskedChangePoints:oldmask");
    integer numchanged = 0;
    for(integer q=1; q<nPoints; q++)
    {
        integer p = SortedIndices[q];
        SafeArray<integer> newmask(Masks, p*nDims,
                "ComputeSortedUnmaskedChangePoints:newmask");
        bool changed = false;
        for(integer i=0; i<nDims; i++)
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
    bool operator()(const integer i, const integer j) const;
};

// Less than operator for KK.Masks, it's just a lexicographical comparison
bool KKSort::operator()(const integer i, const integer j) const
{
    integer nDims = kk->nDims;
    for(integer k=0; k<nDims; k++)
    {
        integer x = kk->Masks[i*nDims+k];
        integer y = kk->Masks[j*nDims+k];
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
    for(integer i=0; i<nPoints; i++)
        SortedIndices[i] = i;
    stable_sort(SortedIndices.begin(), SortedIndices.end(), kksorter);
}

//void KK::ComputeNoiseMeansAndVariances()
//{
    //For TrySplits
    //maintain noise mean and variance of each channel
    // maintain number of masked points in each channel
//    Output("ComputeNoiseMeansandVariances ");

//    NoiseMean.resize(nDims);
//    NoiseVariance.resize(nDims);
//    nMasked.resize(nDims);

//}


void KK::ComputeNoiseMeansAndVariances()
{
    // compute noise mean and variance of each channel
    // compute number of masked points in each channel
    Output("Masked EM: Computing Noise Means and Variances \n -----------------------------------------");

    NoiseMean.resize(nDims);
    NoiseVariance.resize(nDims);
    nMasked.resize(nDims);
    
//    for(integer i=0; i<nDims; i++)
//    {    NoiseMean[i] = 0;
//        NoiseVariance[i] = 0;
//        nMasked[i] = 0;
//    }
        
    
    for(integer p=0; p<nPoints; p++)
        for(integer i=0; i<nDims; i++)
            if(!Masks[p*nDims+i])
            {
                scalar thisdata = Data[p*nDims+i];
                NoiseMean[i] += thisdata;
            //    NoiseVariance[i] += thisdata*thisdata; // sum of squares
                nMasked[i]++;
            }
    for(integer i=0; i<nDims; i++)
    {
        if(nMasked[i]==0)
        {
            NoiseMean[i] = 0.0;
            NoiseVariance[i] = 0;
    //        NoiseVariance[i] = 1.0;
        } else
        {
            NoiseMean[i] /= (scalar)nMasked[i];
        //    NoiseVariance[i] /= (scalar)nMasked[i]; // E[X^2]
        //    NoiseVariance[i] -= NoiseMean[i]*NoiseMean[i]; // -E[X]^2
        }
    }
    
    
    for(integer p=0; p<nPoints; p++)
        for(integer i=0; i<nDims; i++)
            if(!Masks[p*nDims+i])
            {    scalar thisdata = Data[p*nDims+i];
                NoiseVariance[i] += (thisdata-NoiseMean[i])*(thisdata-NoiseMean[i]); 
                
            }
    
    for(integer i=0; i<nDims; i++)
    {   
        if(nMasked[i]==0)
        {    NoiseVariance[i]= 0;
            
        }else {
            NoiseVariance[i] /= (scalar)nMasked[i];
        }

    }
    
//    for(integer i=0; i<nDims; i++)
//    {    Output(" NoiseMean[%d] = %f",i,NoiseMean[i]);
//        Output(" NoiseVariance[%d] = %f",i,NoiseVariance[i]);
//        Output(" nMasked[%d] = %d",i,nMasked[i]);
        
//    }
}

void KK::ComputeCorrectionTermsAndReplaceData()
{
    for(integer p=0; p<nPoints; p++)
        for(integer i=0; i<nDims; i++)
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
    integer i,p;
    
    
    
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
