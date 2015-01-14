/*
 * Main header file
 */

#ifndef MASKED_KLUSTA_KWIK_2_H_
#define MASKED_KLUSTA_KWIK_2_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include "parameters.h"
#include "log.h"
#include "util.h"
#include "numerics.h"
#include "linalg.h"
#include "memorytracking.h"

using namespace std;

class CompoundScore
{
public:
	scalar raw;
	scalar total;
	scalar penalty;
	CompoundScore() : raw(0.0), total(0.0), penalty(0.0) {};
	CompoundScore(scalar _raw, scalar _total, scalar _penalty) : raw(_raw), total(_total), penalty(_penalty) {};
};

class KK {
public:
    /////////////// CONSTRUCTORS ///////////////////////////////////////////////
    // Construct from file
    KK(char *FileBase, integer ElecNo, char *UseFeatures,
            scalar PenaltyK, scalar PenaltyKLogN, integer PriorPoint);
    // Construct from subset of existing KK object
    void ConstructFrom(const KK &Source, const vector<integer> &Indices);
	void ConstructFrom(const KK &Source);
    KK(const KK &Source, const vector<integer> &Indices);
    // Make an entire copy of existing KK object
    KK(const KK &Source);
	// Destructor
	~KK();
    /////////////// FUNCTIONS //////////////////////////////////////////////////
	void MemoryCheck();
    void AllocateArrays();
    void Reindex();
    // Random initial starting conditions functions
    void StartingConditionsRandom();
    void StartingConditionsFromMasks();
    // Precomputation related functions
    void DoInitialPrecomputations();
    // Functions called by DoPrecomputations
    void DoPrecomputations();
    //SNK Precomputation for TrySplits to avoid changing NoiseMean and NoiseVariance and nMasked
    void ComputeUnmasked();
    void ComputeSortIndices();
    void ComputeSortedUnmaskedChangePoints();
    void ComputeNoiseMeansAndVariances();
    void ComputeCorrectionTermsAndReplaceData();
    void PointMaskDimension ();//SNK PointMaskDimension() computes the sum of the masks/float masks for each point
    // Precomputations for cluster masks
    void ComputeClusterMasks();
    // Score and penalty functions
    CompoundScore ComputeScore();
    scalar Penalty(integer n);
    void ComputeClassPenalties();
    // Main algorithm functions
    void MStep();
    void EStep();
    void CStep(bool allow_assign_to_noise=true);
    void ConsiderDeletion();
    integer TrySplits();
	CompoundScore CEM(char *CluFile, integer recurse, integer InitRand, bool allow_assign_to_noise = true);
	CompoundScore Cluster(char *CluFile);
    // IO related functions
    void LoadData(char *FileBase, integer ElecNo, char *UseFeatures);
    void LoadClu(char *StartCluFile);
    void SaveOutput();
    void SaveTempOutput();
    void SaveSortedData();
    void SaveSortedClu();
    void SaveCovMeans();
public:
    /////////////// VARIABLES //////////////////////////////////////////////////
	KK *KK_split, *K2_container; // used for splitting
    integer nDims, nDims2; // nDims2 is nDims squared and the mean of the unmasked dimensions.
    int nStartingClusters; // total # starting clusters, including clu 0, the noise cluster (int because read in from file)
    integer nClustersAlive; // nClustersAlive is total number with points in, excluding noise cluster
    integer nPoints;
    integer priorPoint; // Prior for regularization NOTE: separate from global variabl PriorPoint (capitalization)
    integer NoisePoint; // number of fake points always in noise cluster to ensure noise weight>0
    integer FullStep; // Indicates that the next E-step should be a full step (no time saving)
    scalar penaltyK, penaltyKLogN;

#ifdef STORE_DATA_AS_INTEGER
    vector<data_int> Data; // Data[p*nDims + d] = Input data for point p, dimension d
#else
    vector<scalar> Data; // Data[p*nDims + d] = Input data for point p, dimension d
#endif
    
    // We sort the points into an order where the corresponding mask changes as
    // infrequently as possible, this vector is used to store the sorted indices,
    // see the function KK::ComputeSortIndices() in sortdata.cpp to see how this
    // sorting works. In the case of the TrySplits() function this sorting needs
    // to be recomputed when we change the underlying data.
    vector<integer> SortedIndices;

    vector<char> Masks; //SNK: Masks[p*nDims + d] = Input masks for point p, dimension d
#ifdef COMPUTED_BINARY_MASK
	inline char GetMasks(integer i)
	{
		if(!UseDistributional)
			return Masks[i];
#ifdef STORE_FLOAT_MASK_AS_CHAR
		return CharFloatMasks[i]==(unsigned char)255;
#else
		return FloatMasks[i]==(scalar)1;
#endif
	}
#else
	inline char GetMasks(integer i) { return Masks[i]; }
#endif

#ifdef STORE_FLOAT_MASK_AS_CHAR
	vector<unsigned char> CharFloatMasks; // float mask that is stored in a char to save RAM
#else
    vector<scalar> FloatMasks; // as above but for floating point masks
#endif
    // We store just the indices of the unmasked points in this sparse array
    // structure. For point p, the segment Unmasked[UnmaskedInd[p]] to
    // Unmasked[UnmaskedInd[p+1]] contains the indices i where Masks[i]==1.
    // This allows for a nice contiguous piece of memory without a separate
    // memory allocation for each point (and there can be many points). The
    // precomputation is performed by the function ComputeUnmasked(), and is
    // automatically called by LoadData() at the appropriate time, although
    // in the case of TrySplits() it has to be called explicitly.
    vector<integer> Unmasked;
    vector<integer> UnmaskedInd;
    // We also store a bool that tells us if the mask has changed when we go
    // through the points in sorted order. Computed by function
    // ComputeSortedUnmaskedChangePoints()
    vector<integer> SortedMaskChange;
    
    vector<scalar> UnMaskDims; //SNK: Number of unmasked dimensions for each data point. In the Float masks case, the sum of the weights.
    vector<integer> nClassMembers;
    
    vector<scalar> Weight; // Weight[c] = Class weight for class c
    vector<scalar> Mean; // Mean[c*nDims + d] = cluster mean for cluster c in dimension d
    vector<scalar> Cov; // Cov[c*nDims*nDims + i*nDims + j] = Covariance for cluster C, entry i,j
                    // NB covariances are stored in upper triangle (j>=i)
	vector<BlockPlusDiagonalMatrix> DynamicCov; // Covariance matrices, DynamicCov[cc] where cc is alive cluster index
    vector<scalar> LogP; // LogP[p*MaxClusters + c] = minus log likelihood for point p in cluster c
    vector<integer> Class; // Class[p] = best cluster for point p
    vector<integer> OldClass; // Class[p] = previous cluster for point p
    vector<integer> Class2; // Class[p] = second best cluster for point p
    vector<integer> BestClass; // BestClass = best classification yet achieved
    vector<integer> ClassAlive; // contains 1 if the class is still alive - otherwise 0
    vector<integer> AliveIndex; // a list of the alive classes to iterate over

    // Used for distributional optimisations
    vector<scalar> ClusterMask;
    vector< vector<integer> > ClusterUnmaskedFeatures;
    vector< vector<integer> > ClusterMaskedFeatures;

    // Used in EStep(), but this will probably change later
    vector<scalar> AllVector2Mean;
    // used in M step
    vector<scalar> NoiseMean;
    vector<scalar> NoiseVariance;
    vector<integer> nMasked;
    // used in distribution EM steps

#ifdef COMPUTED_CORRECTION_TERM
	inline scalar GetCorrectionTerm(integer p, integer i)
	{
            // scalar x = Data[p*nDims+i];
#ifdef STORE_FLOAT_MASK_AS_CHAR
            scalar w = CharFloatMasks[p*nDims+i]/(scalar)255.0;
#else
            scalar w = FloatMasks[p*nDims+i];
#endif
            scalar nu = NoiseMean[i];
            scalar sigma2 = NoiseVariance[i];
            //scalar y = w*x+(1-w)*nu;
            //scalar z = w*x*x+(1-w)*(nu*nu+sigma2);
			scalar y = GetData(p, i);
			if(w==(scalar)0.0)
			{
				scalar z = nu*nu+sigma2;
				return z-y*y;
			} else
			{
				scalar x = (y-(1-w)*nu)/w;
				scalar z = w*x*x+(1-w)*(nu*nu+sigma2);
				return z-y*y;
			}
	};
#else
    vector<scalar> CorrectionTerm;
#endif
#ifdef STORE_DATA_AS_INTEGER
	inline scalar GetData(integer p, integer i)
	{
		return data_scalar_from_int(Data[p*nDims+i]);
	};
#else
	inline scalar GetData(integer p, integer i)
	{
		return Data[p*nDims+i];
	};
#endif

    // used in ComputeScore and ConsiderDeletion
    vector<scalar> ClassPenalty;
    // debugging info
    integer numiterations;
	integer init_type;
};

#endif /* MASKED_KLUSTA_KWIK_2_H_ */
