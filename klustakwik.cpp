// MaskedKlustaKwik2.C
//
// Fast clustering using the CEM algorithm with Masks.

// Disable some Visual Studio warnings
#define _CRT_SECURE_NO_WARNINGS

#include "klustakwik.h"
#define _USE_MATH_DEFINES
#include<math.h>

// GLOBAL VARIABLES
FILE *Distfp;
integer global_numiterations = 0;
scalar iteration_metric2 = (scalar)0;
scalar iteration_metric3 = (scalar)0;
clock_t Clock0;
scalar timesofar;

// Does a memory check (should only be called for first instance of KK)
void KK::MemoryCheck()
{
	integer num_bytes_required = 3 * NumBytesRequired();
	scalar memory_required = (num_bytes_required*1.0) / (1024.0*1024.0*1024.0);
	if (memory_required > memory_tracker.limit_gb)
	{
		Error("Running KlustaKwik on this data will use between %.2f and %.2f GB of RAM, and the limit is set at %.2f.\n", (double)(memory_required*2.0 / 3.0), (double)memory_required, (double)memory_tracker.limit_gb);
		Error("Possible candidates are:\n");
		Error("- nPoints = %d\n", (int)nPoints);
		Error("- nDims = %d\n", (int)nDims);
		Error("- MaxPossibleClusters = %d\n", (int)MaxPossibleClusters);
		exit(EXIT_FAILURE);
	}
	Output("This run is expected to use between %.2f and %.2f GB of RAM.\n", (double)(memory_required*2.0 / 3.0), (double)memory_required);
}

integer KK::NumBytesRequired()
{
	// we don't allocate any memory if we have already allocated memory to this
	// (i.e. if we are in TrySplits)
	if (Data.size())
		return 0;
	nDims2 = nDims*nDims;
	// Compute required memory and check if it exceeds the limit set
	integer num_bytes_allocated =
		sizeof(scalar)*nPoints*nDims +               // Data
		sizeof(integer)*nPoints*nDims +              // Masks
		sizeof(scalar)*nPoints*nDims +               // FloatMasks
		sizeof(scalar)*nPoints +                     // UnMaskDims
		sizeof(scalar)*MaxPossibleClusters +         // Weight
		sizeof(scalar)*MaxPossibleClusters*nDims +   // Mean
		(1 - UseDistributional)*sizeof(scalar)*MaxPossibleClusters*nDims2 + // Cov
		sizeof(scalar)*MaxPossibleClusters*nPoints + // LogP
		sizeof(integer)*nPoints +                    // Class
		sizeof(integer)*nPoints +                    // OldClass
		sizeof(integer)*nPoints +                    // Class2
		sizeof(integer)*nPoints +                    // BestClass
		sizeof(integer)*MaxPossibleClusters +        // ClassAlive
		sizeof(integer)*MaxPossibleClusters +        // AliveIndex
		sizeof(scalar)*MaxPossibleClusters +         // ClassPenalty
		sizeof(integer)*MaxPossibleClusters +        // nClassMembers
		sizeof(scalar)*nPoints*nDims +               // AllVector2Mean
		// UseDistributional only
		UseDistributional*sizeof(scalar)*MaxPossibleClusters +  // CorrectionTerm
		sizeof(scalar)*(UseDistributional*MaxPossibleClusters*nDims) + // ClusterMask (vector<scalar>)
		UseDistributional*sizeof(integer)*MaxPossibleClusters*nDims; // ClusterUnmaskedFeatures + ClusterMaskedFeatures

	return num_bytes_allocated;
}

template<class T>
inline void resize_and_fill_with_zeros(vector<T> &x, integer newsize)
{
	if (x.size() == 0)
	{
		x.resize((uinteger)newsize);
		return;
	}
	if (x.size() > (uinteger)newsize)
	{
		fill(x.begin(), x.end(), (T)0);
		x.resize((uinteger)newsize);
	}
	else
	{
		x.resize((uinteger)newsize);
		fill(x.begin(), x.end(), (T)0);
	}
}

// Sets storage for KK class.  Needs to have nDims and nPoints defined
void KK::AllocateArrays() {

    nDims2 = nDims*nDims;
    NoisePoint = 1; // Ensures that the mixture weight for the noise cluster never gets to zero

	integer num_bytes_allocated = NumBytesRequired();
	mem.add(num_bytes_allocated);

    // Set sizes for arrays
	resize_and_fill_with_zeros(Data, nPoints * nDims);
    //SNK
	resize_and_fill_with_zeros(Masks, nPoints * nDims);
	resize_and_fill_with_zeros(FloatMasks, nPoints * nDims);
	resize_and_fill_with_zeros(UnMaskDims, nPoints); //SNK Number of unmasked dimensions for each data point when using float masks $\sum m_i$
	resize_and_fill_with_zeros(Weight, MaxPossibleClusters);
	resize_and_fill_with_zeros(Mean, MaxPossibleClusters*nDims);
	if (!UseDistributional)
		resize_and_fill_with_zeros(Cov, MaxPossibleClusters*nDims2);
	resize_and_fill_with_zeros(LogP, MaxPossibleClusters*nPoints);
	resize_and_fill_with_zeros(Class, nPoints);
	resize_and_fill_with_zeros(OldClass, nPoints);
	resize_and_fill_with_zeros(Class2, nPoints);
	resize_and_fill_with_zeros(BestClass, nPoints);
	resize_and_fill_with_zeros(ClassAlive, MaxPossibleClusters);
	resize_and_fill_with_zeros(AliveIndex, MaxPossibleClusters);
	resize_and_fill_with_zeros(ClassPenalty, MaxPossibleClusters);
	resize_and_fill_with_zeros(nClassMembers, MaxPossibleClusters);
    if(UseDistributional)
    {
		resize_and_fill_with_zeros(CorrectionTerm, nPoints * nDims);
		resize_and_fill_with_zeros(ClusterMask, MaxPossibleClusters*nDims);
    }
}

// recompute index of alive clusters (including 0, the noise cluster)
// should be called after anything that changes ClassAlive
void KK::Reindex()
{
    integer c;

    AliveIndex[0] = 0;
    nClustersAlive=1;
    for(c=1;c<MaxPossibleClusters;c++)
    {
        if (ClassAlive[c])
        {
            AliveIndex[nClustersAlive] = c;
            nClustersAlive++;
        }
    }
}



// Penalty for standard CEM
// Penalty(nAlive) returns the complexity penalty for that many clusters
// bearing in mind that cluster 0 has no free params except p.
scalar KK::Penalty(integer n)
{
    integer nParams;
    if(n==1)
        return 0;
    nParams = (nDims*(nDims+1)/2 + nDims + 1)*(n-1); // each has cov, mean, &p 

    scalar p = penaltyK*(scalar)(nParams) // AIC units (Spurious factor of 2 removed from AIC units on 09.07.13)
            +penaltyKLogN*((scalar)nParams*(scalar)log((scalar)nPoints)/2); // BIC units
    return p;
}


// Penalties for Masked CEM
void KK::ComputeClassPenalties()
{
    if(!((bool)UseDistributional)) // This function must only be called in Use Distributional  mode
    {
  //      Output("Caught in ComputeClassPenalties");
        return;
    }
   // Output("ComputeClassPenalties: Correct if UseDistributional only");
    for(integer c=0; c<MaxPossibleClusters; c++)
        ClassPenalty[c] = (scalar)0;
    // compute sum of nParams for each
    vector<integer> NumberInClass(MaxPossibleClusters);
    for(integer p=0; p<nPoints; p++)
    {
        integer c = Class[p];
        NumberInClass[c]++;
    //    integer n = UnmaskedInd[p+1]-UnmaskedInd[p]; // num unmasked dimensions
        scalar n = UnMaskDims[p];
        scalar nParams = n*(n+1)/2+n+1;
        ClassPenalty[c] += nParams;
    }
    // compute mean nParams for each cluster
    for(integer c=0; c<MaxPossibleClusters; c++)
        if(NumberInClass[c]>0)
            ClassPenalty[c] /= (scalar)NumberInClass[c];
    // compute penalty for each cluster
    for(integer c=0; c<MaxPossibleClusters; c++)
    {
        scalar nParams = ClassPenalty[c];
        ClassPenalty[c] = penaltyK*(scalar)(nParams*2)
                +penaltyKLogN*((scalar)nParams*(scalar)log((scalar)nPoints)/2);
    }
}


// Compute the cluster masks (i.e. the sets of features which are masked/unmasked
// for the whole cluster). Used by M-step and E-step
void KK::ComputeClusterMasks()
{
	Reindex();

	// Initialise cluster mask to 0
	for(integer i=0; i<nDims*MaxPossibleClusters; i++)
		ClusterMask[i] = 0;

	// Compute cluster mask
    for(integer p=0; p<nPoints; p++)
    {
        integer c = Class[p];
		for (integer i = 0; i < nDims; i++)
		{
			ClusterMask[c*nDims + i] += FloatMasks[p*nDims + i];
		}
    }

    // Compute the set of masked/unmasked features for each cluster

	// reset all the subvectors to empty
	ClusterUnmaskedFeatures.clear();
	ClusterUnmaskedFeatures.resize(MaxPossibleClusters);
	ClusterMaskedFeatures.clear();
	ClusterMaskedFeatures.resize(MaxPossibleClusters);
	// fill them in
	for (integer cc = 0; cc<nClustersAlive; cc++)
	{
		integer c = AliveIndex[cc];
		vector<integer> &CurrentUnmasked = ClusterUnmaskedFeatures[c];
		vector<integer> &CurrentMasked = ClusterMaskedFeatures[c];
		for (integer i = 0; i < nDims; i++)
		{
			if (ClusterMask[c*nDims + i]>=PointsForClusterMask)
				CurrentUnmasked.push_back(i);
			else
				CurrentMasked.push_back(i);
		}
		if (Verbose>=2)
		{
			Output("Cluster mask: cluster %d unmasked %d iterations %d/%d init type %d.\n",
				(int)cc, (int)CurrentUnmasked.size(),
				(int)numiterations, (int)global_numiterations, (int)init_type);
		}
	}

}


// M-step: Calculate mean, cov, and weight for each living class
// also deletes any classes with fewer points than nDim
void KK::MStep()
{
    integer p, c, cc, i, j;
    vector<scalar> Vec2Mean(nDims);

    // clear arrays
    memset((void*)&nClassMembers.front(), 0, MaxPossibleClusters*sizeof(integer));
    memset((void*)&Mean.front(), 0, MaxPossibleClusters*nDims*sizeof(scalar));
	if (!UseDistributional)
		memset((void*)&Cov.front(), 0, MaxPossibleClusters*nDims2*sizeof(scalar));
// NOTE: memset commands above replace the code below:
//    for(c=0; c<MaxPossibleClusters; c++) {
//        nClassMembers[c] = 0;
//        for(i=0; i<nDims; i++) Mean[c*nDims + i] = 0;
//    }

    if (Debug) { Output("Entering Unmasked Mstep \n");}

    // Accumulate total number of points in each class
    for (p=0; p<nPoints; p++) nClassMembers[Class[p]]++;

    // check for any dead classes

    if(UseDistributional)
    {
        for (integer cc=0; cc<nClustersAlive; cc++)
        {
            integer c = AliveIndex[cc];
            if (Debug){Output("DistributionalMstep: Class %d contains %d members \n", (int)c, (int)nClassMembers[c]);}
                if (c>0 && nClassMembers[c]<1)//nDims)
                {
                    ClassAlive[c]=0;
                    if (Debug) {Output("UnmaskedMstep_dist: Deleted class %d: no members\n", (int)c);}
                }
        }
    }
    else
    {

        for (cc=0; cc<nClustersAlive; cc++)
        {
                c = AliveIndex[cc];
                if (Debug) {Output("Mstep: Class %d contains %d members \n", (int)c, (int)nClassMembers[c]);}
                if (c>0 && nClassMembers[c]<=nDims)
                {
                    ClassAlive[c]=0;
                    if (Debug) {Output("Deleted class %d: not enough members\n", (int)c);}
                }
        }
    }


    Reindex();


    // Normalize by total number of points to give class weight
    // Also check for dead classes
    if(UseDistributional)
    {
        for (cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
            //Output("DistributionalMstep: PriorPoint on weights ");
            // add "noise point" to make sure Weight for noise cluster never gets to zero
            if(c==0)
            {
                Weight[c] = ((scalar)nClassMembers[c]+NoisePoint) / (nPoints+NoisePoint+priorPoint*(nClustersAlive-1));
            }
            else
            {
                Weight[c] = ((scalar)nClassMembers[c]+priorPoint) / (nPoints+NoisePoint+priorPoint*(nClustersAlive-1));
            }
        }
    }
    else  // For Original KlustaKwik, Classical EM
    {
        for (cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
            // add "noise point" to make sure Weight for noise cluster never gets to zero
            if(c==0)
            {
                Weight[c] = ((scalar)nClassMembers[c]+NoisePoint) / (nPoints+NoisePoint);
            }
            else
            {
                Weight[c] = ((scalar)nClassMembers[c]) / (nPoints+NoisePoint);
            }
        }

    }

    Reindex();

    // Accumulate sums for mean calculation
    for (p=0; p<nPoints; p++)
    {
        c = Class[p];
        for(i=0; i<nDims; i++)
        {
            Mean[c*nDims + i] += Data[p*nDims + i];
        }
    }

    // and normalize
    for (cc=0; cc<nClustersAlive; cc++)
    {
        c = AliveIndex[cc];
        for (i=0; i<nDims; i++) Mean[c*nDims + i] /= nClassMembers[c];
    }

    // Covariance matrix is quite big, and won't fit in the L1d cache
    // (which is 16 or 32 k usually, corresponding to a matrix of about 64x64 or 90x90)
    // so can probably improve performance by doing some sort of blocking here
    // Accumulate sums for covariance calculation
//    for (p=0; p<nPoints; p++)
//    {
//        c = Class[p];
//        // calculate distance from mean
//        for(i=0; i<nDims; i++)
//            Vec2Mean[i] = Data[p*nDims + i] - Mean[c*nDims + i];
//        for(i=0; i<nDims; i++)
//            for(j=i; j<nDims; j++)
//                Cov[c*nDims2 + i*nDims + j] += Vec2Mean[i] * Vec2Mean[j];
//    }
	if ((integer)AllVector2Mean.size() < nPoints*nDims)
	{
		//mem.add((nPoints*nDims-AllVector2Mean.size())*sizeof(scalar));
		AllVector2Mean.resize(nPoints*nDims);
	}
    vector< vector<integer> > PointsInClass(MaxPossibleClusters);
    for(p=0; p<nPoints; p++)
    {
        c = Class[p];
        PointsInClass[c].push_back(p);
		for (i = 0; i < nDims; i++)
			AllVector2Mean[p*nDims + i] = Data[p*nDims + i] - Mean[c*nDims + i];
    }

	if (UseDistributional)
	{
		// Compute the cluster masks, used below to optimise the computation
		ComputeClusterMasks();
		// Empty the dynamic covariance matrices (we will fill it up as we go)
		DynamicCov.clear();

		for (cc = 0; cc<nClustersAlive; cc++)
		{
			c = AliveIndex[cc];
			vector<integer> &PointsInThisClass = PointsInClass[c];
			vector<integer> &CurrentUnmasked = ClusterUnmaskedFeatures[c];
			vector<integer> &CurrentMasked = ClusterMaskedFeatures[c];
			DynamicCov.push_back(BlockPlusDiagonalMatrix(CurrentMasked, CurrentUnmasked));
			BlockPlusDiagonalMatrix &CurrentCov = DynamicCov.back();
			if (CurrentUnmasked.size() == 0)
				continue;

			//// Correct version for dynamic cov matrix
			//for (integer q = 0; q < (integer)PointsInThisClass.size(); q++)
			//{
			//	p = PointsInThisClass[q];
			//	for (integer ii = 0; ii < (integer)CurrentUnmasked.size(); ii++)
			//	{
			//		i = CurrentUnmasked[ii];
			//		for (integer jj = 0; jj < (integer)CurrentUnmasked.size(); jj++)
			//		{
			//			j = CurrentUnmasked[jj];
			//			//Cov[c*nDims2 + i*nDims + j] += AllVector2Mean[p*nDims + i] * AllVector2Mean[p*nDims + j];
			//			CurrentCov.Block[ii*CurrentCov.NumUnmasked + jj] += AllVector2Mean[p*nDims + i] * AllVector2Mean[p*nDims + j];
			//		}
			//	}
			//}
			// Fast version for dynamic cov matrix
			const integer npoints = (integer)PointsInThisClass.size();
			const integer nunmasked = (integer)CurrentUnmasked.size();
			if (npoints > 0 && nunmasked > 0)
			{
				const integer * __restrict pitc = &(PointsInThisClass[0]);
				const integer * __restrict cu = &(CurrentUnmasked[0]);
				for (integer q = 0; q < npoints; q++)
				{
					const integer p = pitc[q];
					const scalar * __restrict av2mp = &(AllVector2Mean[p*nDims]);
					for (integer ii = 0; ii < nunmasked; ii++)
					{
						const integer i = cu[ii];
						const scalar av2mp_i = av2mp[i];
						scalar * __restrict row = &(CurrentCov.Block[ii*nunmasked]);
						for (integer jj = 0; jj < nunmasked; jj++)
						{
							const integer j = cu[jj];
							//Cov[c*nDims2 + i*nDims + j] += AllVector2Mean[p*nDims + i] * AllVector2Mean[p*nDims + j];
							row[jj] += av2mp_i * av2mp[j];
						}
					}
				}
			}

			// Correct version
			//for (integer q = 0; q < (integer)PointsInThisClass.size(); q++)
			//{
			//	p = PointsInThisClass[q];
			//	for (integer ii = 0; ii < (integer)CurrentUnmasked.size(); ii++)
			//	{
			//		i = CurrentUnmasked[ii];
			//		for (integer jj = 0; jj < (integer)CurrentUnmasked.size(); jj++)
			//		{
			//			j = CurrentUnmasked[jj];
			//			Cov[c*nDims2 + i*nDims + j] += AllVector2Mean[p*nDims + i] * AllVector2Mean[p*nDims + j];
			//		}
			//	}
			//}
			// Faster version (equivalent)
			// Doesn't make any use of cache structure, but no need to upgrade now because
			// we will move to a sparse block matrix structure that will make this more
			// natural
			/*
			const integer * __restrict cu = &(CurrentUnmasked[0]);
			const integer ncu = (integer)CurrentUnmasked.size();
			scalar * __restrict cov_c = &(Cov[c*nDims2]);
			const integer * __restrict pitc = &(PointsInThisClass[0]);
			const integer npitc = (integer)PointsInThisClass.size();
			const scalar * __restrict av2m = &(AllVector2Mean[0]);
			for (integer q = 0; q < npitc; q++)
			{
				const integer p = pitc[q];
				const scalar * __restrict av2m_p = av2m + p*nDims;
				for (integer ii = 0; ii < ncu; ii++)
				{
					const integer i = cu[ii];
					const scalar av2m_p_i = av2m_p[i];
					scalar * __restrict cov_c_i = cov_c + i*nDims;
					for (integer jj = 0; jj < ncu; jj++)
					{
						const integer j = cu[jj];
						cov_c_i[j] += av2m_p_i*av2m_p[j];
						//Cov[c*nDims2 + i*nDims + j] += AllVector2Mean[p*nDims + i] * AllVector2Mean[p*nDims + j];
					}
				}
			}
			*/
		}
	}
	else
	{
		// I think this code gives wrong results (but only slightly) (DFMG: 2014/10/13)
		for (c = 0; c < MaxPossibleClusters; c++)
		{
			vector<integer> &PointsInThisClass = PointsInClass[c];
			SafeArray<scalar> safeCov(Cov, c*nDims2, "safeCovMStep");
			for (integer iblock = 0; iblock < nDims; iblock += COVARIANCE_BLOCKSIZE)
			{
				for (integer jblock = iblock; jblock < nDims; jblock += COVARIANCE_BLOCKSIZE)
				{
					for (integer q = 0; q < (integer)PointsInThisClass.size(); q++)
					{
						p = PointsInThisClass[q];
						scalar *cv2m = &AllVector2Mean[p*nDims];
						for (i = iblock; i < MIN(nDims, iblock + COVARIANCE_BLOCKSIZE); i++)
						{
							scalar cv2mi = cv2m[i];
							integer jstart;
							if (jblock != iblock)
								jstart = jblock;
							else
								jstart = i;
							scalar *covptr = &safeCov[i*nDims + jstart];
							scalar *cv2mjptr = &cv2m[jstart];
							//scalar *cv2mjend = cv2m+MIN(nDims, jblock+COVARIANCE_BLOCKSIZE);
							//for(j=jstart; j<MIN(nDims, jblock+COVARIANCE_BLOCKSIZE); j++)
							//for(; cv2mjptr!=cv2mjend;)
							for (j = MIN(nDims, jblock + COVARIANCE_BLOCKSIZE) - jstart; j; j--)
								*covptr++ += cv2mi*(*cv2mjptr++);
						}
					}
				}
			}
		}
	}

    if(UseDistributional)
    {
        for(cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
            vector<integer> &PointsInThisClass = PointsInClass[c];
			integer NumPointsInThisClass = PointsInThisClass.size();
			BlockPlusDiagonalMatrix &CurrentCov = DynamicCov[cc];
			for (integer ii = 0; ii<CurrentCov.NumUnmasked; ii++)
			{
				i = (*CurrentCov.Unmasked)[ii];
				scalar ccf = 0.0; // class correction factor
				for (integer q = 0; q<NumPointsInThisClass; q++)
				{
					p = PointsInThisClass[q];
					ccf += CorrectionTerm[p*nDims + i];
				}
				CurrentCov.Block[ii*CurrentCov.NumUnmasked + ii] += ccf;
			}
			for (integer ii = 0; ii<CurrentCov.NumMasked; ii++)
			{
				i = (*CurrentCov.Masked)[ii];
				scalar ccf = 0.0; // class correction factor
				for (integer q = 0; q<NumPointsInThisClass; q++)
				{
					p = PointsInThisClass[q];
					ccf += CorrectionTerm[p*nDims + i];
				}
				CurrentCov.Diagonal[ii] += ccf;
			}
		}

    // Add a diagonal matrix of Noise variances to the covariance matrix for renormalization
		for (cc = 0; cc < nClustersAlive; cc++)
		{
			c = AliveIndex[cc];
			BlockPlusDiagonalMatrix &CurrentCov = DynamicCov[cc];
			for (integer ii = 0; ii < CurrentCov.NumUnmasked; ii++)
				CurrentCov.Block[ii*CurrentCov.NumUnmasked + ii] += priorPoint*NoiseVariance[(*CurrentCov.Unmasked)[ii]];
			for (integer ii = 0; ii < CurrentCov.NumMasked; ii++)
				CurrentCov.Diagonal[ii] += priorPoint*NoiseVariance[(*CurrentCov.Masked)[ii]];
		}


    }

    // and normalize
    if(UseDistributional)
	{
        for (cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
			BlockPlusDiagonalMatrix &CurrentCov = DynamicCov[cc];
			scalar factor = 1.0 / (nClassMembers[c] + priorPoint - 1);
			for (i = 0; i < (integer)CurrentCov.Block.size(); i++)
				CurrentCov.Block[i] *= factor;
			for (i = 0; i < (integer)CurrentCov.Diagonal.size(); i++)
				CurrentCov.Diagonal[i] *= factor;
		}

    }
    else
	{    //For original KlustaKwik classical EM
        for (cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
            for(i=0; i<nDims; i++)
                for(j=i; j<nDims; j++)
                    Cov[c*nDims2 + i*nDims + j] /= (nClassMembers[c]-1);
        }
    }


    // That's it!

    // Diagnostics
    if (Debug)
    {
        for (cc=0; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
			Output("Class %d - Weight %.2g\n", (int)c, Weight[c]);
            Output("Mean: ");
            MatPrint(stdout, &Mean.front() + c*nDims, 1, nDims);
			if (!UseDistributional)
			{
				Output("\nCov:\n");
				MatPrint(stdout, &Cov.front() + c*nDims2, nDims, nDims);
			}
            Output("\n");
        }
    }
}



// E-step.  Calculate Log Probs for each point to belong to each living class
// will delete a class if covariance matrix is singular
// also counts number of living classes
void KK::EStep()
{
    integer p, c, cc, i;
    integer nSkipped;
    scalar LogRootDet; // log of square root of covariance determinant
    scalar Mahal; // Mahalanobis distance of point from cluster center
    scalar correction_factor = (scalar)1; // for partial correction in distributional step
	scalar InverseClusterNorm;
	vector<scalar> Chol(nDims2); // to store choleski decomposition
    vector<scalar> Vec2Mean(nDims); // stores data point minus class mean
    vector<scalar> Root(nDims); // stores result of Chol*Root = Vec
    vector<scalar> InvCovDiag;
    if(UseDistributional)
        InvCovDiag.resize(nDims);

    SafeArray<scalar> safeChol(Chol, "safeChol");
    SafeArray<scalar> safeVec2Mean(Vec2Mean, "safeVec2Mean");
    SafeArray<scalar> safeRoot(Root, "safeRoot");
    SafeArray<scalar> safeInvCovDiag(InvCovDiag, "safeInvCovDiag");

    nSkipped = 0;

    if (Debug) {Output("Entering Unmasked Estep \n");}

    // start with cluster 0 - uniform distribution over space
    // because we have normalized all dims to 0...1, density will be 1.
    vector<integer> NumberInClass(MaxPossibleClusters);  // For finding number of points in each class
    for (p=0; p<nPoints; p++)
    {
        LogP[p*MaxPossibleClusters + 0] = (float)-log(Weight[0]);
        integer ccc = Class[p];
        NumberInClass[ccc]++;
    }

	BlockPlusDiagonalMatrix *CurrentCov;
	BlockPlusDiagonalMatrix *CholBPD = NULL;

	for (cc = 1; cc<nClustersAlive; cc++)
    {
        c = AliveIndex[cc];

        // calculate cholesky decomposition for class c
		integer chol_return;
		if (UseDistributional)
		{
			CurrentCov = &(DynamicCov[cc]);
			if (CholBPD)
			{
				delete CholBPD;
				CholBPD = NULL;
			}
			CholBPD = new BlockPlusDiagonalMatrix(*(CurrentCov->Masked), *(CurrentCov->Unmasked));
			chol_return = BPDCholesky(*CurrentCov, *CholBPD);
			//if (MinMaskOverlap>0)
			//{
			//	// compute the norm of the cluster mask (used for skipping points)
			//	const scalar * __restrict cm = &(ClusterMask[c*nDims]);
			//	scalar ClusterNorm = 0.0;
			//	for (i = 0; i < nDims; i++)
			//	{
			//		scalar m = cm[i];
			//		//if (m > ClusterNorm)
			//		//	ClusterNorm = m;
			//		ClusterNorm += m*m;
			//	}
			//	//InverseClusterNorm = 1.0 / ClusterNorm;
			//	InverseClusterNorm = 1.0 / sqrt(ClusterNorm);
			//	//InverseClusterNorm = sqrt((scalar)nDims) / sqrt(ClusterNorm);
			//}
		}
		else
		{
			SafeArray<scalar> safeCov(Cov, c*nDims2, "safeCov");
			chol_return = Cholesky(safeCov, safeChol, nDims);
		}
        if(chol_return)
        {
            // If Cholesky returns 1, it means the matrix is not positive definite.
            // So kill the class.
            // Cholesky is defined in linalg.cpp
			Output("Unmasked E-step: Deleting class %d (%d points): covariance matrix is singular \n", (int)c, (int)NumberInClass[c]);
            ClassAlive[c] = 0;
            continue;
        }

        // LogRootDet is given by log of product of diagonal elements
		if (UseDistributional)
		{
			LogRootDet = 0;
			for (integer ii = 0; ii < CholBPD->NumUnmasked; ii++)
				LogRootDet += log(CholBPD->Block[ii*CholBPD->NumUnmasked + ii]);
			for (integer ii = 0; ii < CholBPD->NumMasked; ii++)
				LogRootDet += log(CholBPD->Diagonal[ii]);
		}
		else
		{
			LogRootDet = 0;
			for (i = 0; i < nDims; i++)
				LogRootDet += log(Chol[i*nDims + i]);
		}

        // if distributional E step, compute diagonal of inverse of cov matrix
        if(UseDistributional)
        {
            vector<scalar> BasisVector(nDims);
            SafeArray<scalar> safeBasisVector(BasisVector, "BasisVector");
            for(integer i=0; i<nDims; i++)
                safeBasisVector[i] = (scalar)0;
            for(integer i=0; i<nDims; i++)
            {  
				safeBasisVector[i] = (scalar)1;
                // calculate Root vector - by Chol*Root = BasisVector
				BPDTriSolve(*CholBPD, safeBasisVector, safeRoot);
                // add half of Root vector squared to log p
                scalar Sii = (scalar)0;
                for(integer j=0; j<nDims; j++)
                    Sii += Root[j]*Root[j];
                safeInvCovDiag[i] = Sii;
                safeBasisVector[i] = (scalar)0;
            }
        }

        for(p=0; p<nPoints; p++)
        {
            // to save time -- only recalculate if the last one was close
            if (
                !FullStep
                && (Class[p] == OldClass[p])
                && (LogP[p*MaxPossibleClusters+c] - LogP[p*MaxPossibleClusters+Class[p]] > DistThresh)
                )
            {
                nSkipped++;
                continue;
            }

			// to save time, skip points with mask overlap below threshold
			if (MinMaskOverlap > 0)
			{
				// compute dot product of point mask with cluster mask
				const scalar * __restrict PointMask = &(FloatMasks[p*nDims]);
				const scalar * __restrict cm = &(ClusterMask[c*nDims]);
				scalar dotprod = 0.0;
				//// InverseClusterNorm is computed above, uncomment it if you uncomment any of this
				//for (i = 0; i < nDims; i++)
				//{
				//	dotprod += cm[i] * PointMask[i] * InverseClusterNorm;
				//	if (dotprod >= MinMaskOverlap)
				//		break;
				//}
				const integer NumUnmasked = CurrentCov->NumUnmasked;
				if (NumUnmasked)
				{
					const integer * __restrict cu = &((*(CurrentCov->Unmasked))[0]);
					for (integer ii = 0; ii < NumUnmasked; ii++)
					{
						const integer i = cu[ii];
						dotprod += PointMask[i];
						if (dotprod >= MinMaskOverlap)
							break;
					}
				}
				//dotprod *= InverseClusterNorm;
				if (dotprod < MinMaskOverlap)
				{
					nSkipped++;
					continue;
				}
			}

            // Compute Mahalanobis distance
            Mahal = 0;

			// calculate data minus class mean
			//for (i = 0; i<nDims; i++)
			//	Vec2Mean[i] = Data[p*nDims + i] - Mean[c*nDims + i];
			scalar * __restrict Data_p = &(Data[p*nDims]);
			scalar * __restrict Mean_c = &(Mean[c*nDims]);
			scalar * __restrict v2m = &(Vec2Mean[0]);
			for (i = 0; i < nDims; i++)
				v2m[i] = Data_p[i] - Mean_c[i];

            // calculate Root vector - by Chol*Root = Vec2Mean
			if (UseDistributional)
				BPDTriSolve(*CholBPD, safeVec2Mean, safeRoot);
			else
				TriSolve(safeChol, safeVec2Mean, safeRoot, nDims);

            // add half of Root vector squared to log p
            for(i=0; i<nDims; i++)
                Mahal += Root[i]*Root[i];

            // if distributional E step, add correction term
			if (UseDistributional)
			{
				scalar * __restrict ctp = &(CorrectionTerm[p*nDims]);
				scalar * __restrict icd = &(InvCovDiag[0]);
				scalar subMahal = 0.0;
				for (i = 0; i < nDims; i++)
					subMahal += ctp[i] * icd[i];
				Mahal += subMahal*correction_factor;
			}
            // Score is given by Mahal/2 + log RootDet - log weight
            LogP[p*MaxPossibleClusters + c] = Mahal/2
                                       + LogRootDet
                                    - log(Weight[c])
									+ (0.5*log(2 * M_PI))*nDims;

        } // for(p=0; p<nPoints; p++)
    } // for(cc=1; cc<nClustersAlive; cc++)
	if (CholBPD)
		delete CholBPD;
}



// Choose best class for each point (and second best) out of those living
void KK::CStep(bool allow_assign_to_noise)
{
    integer p, c, cc, TopClass, SecondClass;
    integer ccstart = 0;
    if(!allow_assign_to_noise)
        ccstart = 1;
    scalar ThisScore, BestScore, SecondScore;

    for (p=0; p<nPoints; p++)
    {
        OldClass[p] = Class[p];
        BestScore = HugeScore;
        SecondScore = HugeScore;
        TopClass = SecondClass = 0;
        for (cc=ccstart; cc<nClustersAlive; cc++)
        {
            c = AliveIndex[cc];
            ThisScore = LogP[p*MaxPossibleClusters + c];
            if (ThisScore < BestScore)
            {
                SecondClass = TopClass;
                TopClass = c;
                SecondScore = BestScore;
                BestScore = ThisScore;
            }
            else if (ThisScore < SecondScore)
            {
                SecondClass = c;
                SecondScore = ThisScore;
            }
        }
        Class[p] = TopClass;
        Class2[p] = SecondClass;
    }
}

// Sometimes deleting a cluster will improve the score, when you take into account
// the BIC. This function sees if this is the case.  It will not delete more than
// one cluster at a time.
void KK::ConsiderDeletion()
{

    integer c, p, CandidateClass=0;
    scalar Loss, DeltaPen;
    vector<scalar> DeletionLoss(MaxPossibleClusters); // the increase in log P by deleting the cluster
    
    if (Debug)
        Output(" Entering ConsiderDeletion: ");

    for(c=1; c<MaxPossibleClusters; c++)
    {
        if (ClassAlive[c]) DeletionLoss[c] = 0;
        else DeletionLoss[c] = HugeScore; // don't delete classes that are already there
    }

    // compute losses by deleting clusters
    vector<integer> NumberInClass(MaxPossibleClusters);
    for(p=0; p<nPoints; p++)
    {
        DeletionLoss[Class[p]] += LogP[p*MaxPossibleClusters + Class2[p]] - LogP[p*MaxPossibleClusters + Class[p]];
        integer ccc = Class[p];
        NumberInClass[ccc]++;  // For computing number of points in each class
    }

    // find class with smallest increase in total score
    Loss = HugeScore;
    if (UseDistributional) //For UseDistribution, we use the ClusterPenalty
    {
        for(c=1; c<MaxPossibleClusters; c++)
        {
            if ((DeletionLoss[c]-ClassPenalty[c])<Loss)
            {
                Loss = DeletionLoss[c]-ClassPenalty[c];
                CandidateClass = c;
            }
        }
        
    }// or in the case of fixed penalty find class with least to lose
    else
    {
        for(c=1; c<MaxPossibleClusters; c++)
        {
            if (DeletionLoss[c]<Loss)
            {
                Loss = DeletionLoss[c];
                CandidateClass = c;
            }
        }
        
    }

    

    // what is the change in penalty?
    if(UseDistributional) //For the distributional algorithm we need to use the ClusterPenalty
        DeltaPen = ClassPenalty[CandidateClass];
    else
        DeltaPen = Penalty(nClustersAlive) - Penalty(nClustersAlive-1);

    //Output("cand Class %d would lose " SCALARFMT " gain is " SCALARFMT "\n", (int)CandidateClass, Loss, DeltaPen);
    // is it worth it?
    //06/12/12 fixing bug introduced which considered DeltaPen twice!
    if (UseDistributional) //For the distributional algorithm we need to use the ClusterPenalty
    {
        if (Loss<0)
        {
			Output("Deleting Class %d (%d points): Lose " SCALARFMT " but Gain " SCALARFMT "\n", (int)CandidateClass, (int)NumberInClass[CandidateClass], DeletionLoss[CandidateClass], DeltaPen);
            // set it to dead
            ClassAlive[CandidateClass] = 0;
            
            // re-allocate all of its points
            for(p=0;p<nPoints; p++) if(Class[p]==CandidateClass) Class[p] = Class2[p];
            // recompute class penalties
            ComputeClassPenalties();
        }
    }
    else
    {
        if (Loss<DeltaPen)
        {
			Output("Deleting Class %d (%d points): Lose " SCALARFMT " but Gain " SCALARFMT "\n", (int)CandidateClass, (int)NumberInClass[CandidateClass], DeletionLoss[CandidateClass], DeltaPen);
            // set it to dead
            ClassAlive[CandidateClass] = 0;

            // re-allocate all of its points
            for(p=0;p<nPoints; p++) if(Class[p]==CandidateClass) Class[p] = Class2[p];
            // recompute class penalties
            ComputeClassPenalties();
        }
    }

    Reindex();
}


// LoadClu(CluFile)
void KK::LoadClu(char *CluFile)
{
    FILE *fp;
    integer p, c;
	int val; // read in from %d
    integer status;


    fp = fopen_safe(CluFile, "r");
    status = fscanf(fp, "%d", &nStartingClusters);
    nClustersAlive = nStartingClusters;// -1;
    for(c=0; c<MaxPossibleClusters; c++) ClassAlive[c]=(c<nStartingClusters);

    for(p=0; p<nPoints; p++)
    {
        status = fscanf(fp, "%d", &val);
        if (status==EOF) Error("Error reading cluster file");
        Class[p] = val-1;
    }
}

// for each cluster, try to split it in two.  if that improves the score, do it.
// returns 1 if split was successful
integer KK::TrySplits()
{
    integer c, cc, c2, p, p2, DidSplit = 0;
    scalar Score, NewScore, UnsplitScore, SplitScore;
    integer UnusedCluster;
    //KK K2; // second KK structure for sub-clustering
    //KK K3; // third one for comparison

    if(nClustersAlive>=MaxPossibleClusters-1)
    {
        Output("Won't try splitting - already at maximum number of clusters\n");
        return 0;
    }

    // set up K3 and remember to add the masks
    //KK K3(*this);
	if (!AlwaysSplitBimodal)
	{
		if (KK_split == NULL)
		{
			KK_split = new KK(*this);
		}
		else
		{
			// We have to clear these to bypass the debugging checks
			// in precomputations.cpp
			KK_split->Unmasked.clear();
			KK_split->UnmaskedInd.clear();
			KK_split->SortedMaskChange.clear();
			KK_split->SortedIndices.clear();
			// now we treat it as empty
			KK_split->ConstructFrom(*this);
		}
	}
	//KK &K3 = *KK_split;
#define K3 (*KK_split)

    Output("Compute initial score before splitting: ");
    Score = ComputeScore();

    // loop thu clusters, trying to split
    for (cc=1; cc<nClustersAlive; cc++)
    {
        c = AliveIndex[cc];

        // set up K2 structure to contain points of this cluster only

        vector<integer> SubsetIndices;
        for(p=0; p<nPoints; p++)
            if(Class[p]==c)
                SubsetIndices.push_back(p);
        if(SubsetIndices.size()==0)
            continue;

		if (K2_container)
		{
			// We have to clear these to bypass the debugging checks
			// in precomputations.cpp
			K2_container->Unmasked.clear();
			K2_container->UnmaskedInd.clear();
			K2_container->SortedMaskChange.clear();
			K2_container->SortedIndices.clear();
			//K2_container->AllVector2Mean.clear();
			// now we treat it as empty
			K2_container->ConstructFrom(*this, SubsetIndices);
		}
		else
		{
			K2_container = new KK(*this, SubsetIndices);
		}
        //KK K2(*this, SubsetIndices);
		KK &K2 = *K2_container;

        // find an unused cluster
        UnusedCluster = -1;
        for(c2=1; c2<MaxPossibleClusters; c2++)
        {
             if (!ClassAlive[c2])
             {
                 UnusedCluster = c2;
                 break;
             }
        }
        if (UnusedCluster==-1)
        {
            Output("No free clusters, abandoning split");
            return DidSplit;
        }

        // do it
		if (Verbose >= 1) Output("\n Trying to split cluster %d (%d points) \n", (int)c, (int)K2.nPoints);
        K2.nStartingClusters=2; // (2 = 1 clusters + 1 unused noise cluster)
        UnsplitScore = K2.CEM(NULL, 0, 1, false);
        K2.nStartingClusters=3; // (3 = 2 clusters + 1 unused noise cluster)
        SplitScore = K2.CEM(NULL, 0, 1, false);

        // Fix by Michaël Zugaro: replace next line with following two lines
        // if(SplitScore<UnsplitScore) {
        if(K2.nClustersAlive<2) Output("\n Split failed - leaving alone\n");
        if((SplitScore<UnsplitScore)&&(K2.nClustersAlive>=2)) {
			if (AlwaysSplitBimodal)
			{
				DidSplit = 1;
				Output("\n We are always splitting bimodal clusters so it's getting split into cluster %d.\n", (int)UnusedCluster);
				p2 = 0;
				for (p = 0; p < nPoints; p++)
				{
					if (Class[p] == c)
					{
						if (K2.Class[p2] == 1) Class[p] = c;
						else if (K2.Class[p2] == 2) Class[p] = UnusedCluster;
						else Error("split should only produce 2 clusters\n");
						p2++;
					}
					ClassAlive[Class[p]] = 1;
				}
			}
			else
			{
				// will splitting improve the score in the whole data set?

				// assign clusters to K3
				for (c2 = 0; c2 < MaxPossibleClusters; c2++) K3.ClassAlive[c2] = 0;
				//   Output("%d Points in class %d in KKobject K3 ", (int)c2, (int)K3.nClassMembers[c2]);
				p2 = 0;
				for (p = 0; p < nPoints; p++)
				{
					if (Class[p] == c)
					{
						if (K2.Class[p2] == 1) K3.Class[p] = c;
						else if (K2.Class[p2] == 2) K3.Class[p] = UnusedCluster;
						else Error("split should only produce 2 clusters\n");
						p2++;
					}
					else K3.Class[p] = Class[p];
					K3.ClassAlive[K3.Class[p]] = 1;
				}
				K3.Reindex();

				// compute scores

				K3.MStep();
				K3.EStep();
				//Output("About to compute K3 class penalties");
				if (UseDistributional) K3.ComputeClassPenalties(); //SNK Fixed bug: Need to compute the cluster penalty properly, cluster penalty is only used in UseDistributional mode
				NewScore = K3.ComputeScore();
				Output("\nSplitting cluster %d changes total score from " SCALARFMT " to " SCALARFMT "\n", (int)c, Score, NewScore);

				if (NewScore < Score)
				{
					DidSplit = 1;
					Output("\n So it's getting split into cluster %d.\n", (int)UnusedCluster);
					// so put clusters from K3 back into main KK struct (K1)
					for (c2 = 0; c2 < MaxPossibleClusters; c2++) ClassAlive[c2] = K3.ClassAlive[c2];
					for (p = 0; p < nPoints; p++) Class[p] = K3.Class[p];
				}
				else
				{
					Output("\n So it's not getting split.\n");
				}
			}
        }
    }
	return DidSplit;
#undef K3
}

// ComputeScore() - computes total score.  Requires M, E, and C steps to have been run
scalar KK::ComputeScore()
{
    integer p;
   // integer debugadd;

    scalar penalty = (scalar)0;
    if(UseDistributional)  // For distributional algorithm we require the cluster penalty
        for(integer c=0; c<MaxPossibleClusters; c++)
            penalty += ClassPenalty[c];
    else
        penalty = Penalty(nClustersAlive);
    scalar Score = penalty;
    for(p=0; p<nPoints; p++)
    {    //debugadd = LogP[p*MaxPossibleClusters + Class[p]];
        Score += LogP[p*MaxPossibleClusters + Class[p]];
        // Output("point %d: cumulative score " SCALARFMT " adding" SCALARFMT "\n", (int)p, Score, debugadd);
    }
    //Error("Score: " SCALARFMT " Penalty: " SCALARFMT "\n", Score, penalty);
    Output("  Score: Raw " SCALARFMT " + Penalty " SCALARFMT " = " SCALARFMT "\n", Score-penalty, penalty, Score);

    if (Debug) {
        integer c, cc;
        scalar tScore;
        for(cc=0; cc<nClustersAlive; cc++) {
            c = AliveIndex[cc];
            tScore = 0;
            for(p=0; p<nPoints; p++) if(Class[p]==c) tScore += LogP[p*MaxPossibleClusters + Class[p]];
            Output("class %d has subscore " SCALARFMT "\n", c, tScore);
        }
    }

    return Score;
}

// Initialise starting conditions randomly
void KK::StartingConditionsRandom()
{
    // initialize data to random
    if(nStartingClusters>1)
        for(integer p=0; p<nPoints; p++) // No points are put in the noise cluster to begin 
            Class[p] = irand(1, nStartingClusters-1);
    else
        for(integer p=0; p<nPoints; p++) //If there is only one cluster, put all the points in the noise cluster
            Class[p] = 0;

    for(integer c=0; c<MaxPossibleClusters; c++)
        ClassAlive[c] = (c<nStartingClusters);

	if (SplitInfo == 1) Output("\tSP: Assigned %d initial classes randomly.\n", (int)nStartingClusters);
}

// Initialise starting conditions by selecting unique masks at random
void KK::StartingConditionsFromMasks()
{
    integer nClusters2start=0; //SNK To replace nStartingClusters within this variable only
    
    //if (Debug)
    //    Output("StartingConditionsFromMasks: ");
    Output("Starting initial clusters from distinct float masks \n ");
    
    if(nStartingClusters<=1) // If only 1 starting clutser has been requested, assign all the points to cluster 0
    {
        for(integer p=0; p<nPoints; p++)
            Class[p] = 0;
    }
    else   
    {
        integer num_masks = 0;
        for(integer p=0; p<nPoints; p++)
            num_masks += (integer)SortedMaskChange[p];
        
        if((nStartingClusters-1)>num_masks)
        {
            Error("Not enough masks (%d) to generate starting clusters (%d), "
				"so starting with (%d) clusters instead.\n", (int)num_masks,
				(int)nStartingClusters, (int)(num_masks + 1));
            nClusters2start = num_masks+1;
            //return;
        }
        else
        {
            nClusters2start = nStartingClusters;
        }
        // Construct the set of all masks
        vector<bool> MaskUsed;
        vector<integer> MaskIndex(nPoints);
        vector<integer> MaskPointIndex;
        integer current_mask_index = -1;
        for(integer q=0; q<nPoints; q++)
        {
            integer p = SortedIndices[q];
            if(q==0 || SortedMaskChange[p])
            {
                current_mask_index++;
                MaskUsed.push_back(false);
                MaskPointIndex.push_back(p);
            }
            MaskIndex[p] = current_mask_index;
        }
        // Select points at random until we have enough masks
        integer masks_found = 0;
        vector<integer> MaskIndexToUse;
        vector<integer> FoundMaskIndex(num_masks);
        while(masks_found<nClusters2start-1)
        {
            integer p = irand(0, nPoints-1);
            integer mask_index = MaskIndex[p];
            if(!MaskUsed[mask_index])
            {
                MaskIndexToUse.push_back(mask_index);
                MaskUsed[mask_index] = true;
                FoundMaskIndex[mask_index] = masks_found;
                masks_found++;
            }
        }
        // Assign points to clusters based on masks
        for(integer p=0; p<nPoints; p++)
        {
            if(MaskUsed[MaskIndex[p]]) // we included this points mask
                Class[p] = FoundMaskIndex[MaskIndex[p]]+1; // so assign class to mask index
            else // this points mask not included
            {
                // so find closest match
                integer closest_index = 0;
                integer distance = nDims+1;
                vector<integer> possibilities;
                for(integer mi=0; mi<nClusters2start-1; mi++)
                {
                    integer mip = MaskPointIndex[MaskIndexToUse[mi]];
                    // compute mask distance
                    integer curdistance = 0;
                    for(integer i=0; i<nDims; i++)
                        if(Masks[p*nDims+i]!=Masks[mip*nDims+i])
                            curdistance++;
                    if(curdistance<distance)
                    {
                        possibilities.clear();
                        distance = curdistance;
                    }
                    if(curdistance==distance)
                        possibilities.push_back(mi);
                }
                if((MaskStarts > 0) ||AssignToFirstClosestMask)
                    closest_index = possibilities[0];
                else
                    closest_index = possibilities[irand(0, possibilities.size()-1)];
                Class[p] = closest_index+1;
            }
        }
        // print some info
        Output("Assigned %d initial classes from %d unique masks.\n",
			   (int)nClusters2start, (int)num_masks);
        // Dump initial random classes to a file - knowledge of maskstart configuration may be useful
        // TODO: remove this for final version - SNK: actually it is a nice idea to keep this
        char fname[STRLEN];
        FILE *fp;
		sprintf(fname, "%s.initialclusters.%d.clu.%d", FileBase, (int)nClusters2start, (int)ElecNo);
        fp = fopen_safe(fname, "w");
		fprintf(fp, "%d\n", (int)nClusters2start);
        for(integer p=0; p<nPoints; p++)
			fprintf(fp, "%d\n", (int)Class[p]);
        fclose(fp);
    }
    for(integer c=0; c<MaxPossibleClusters; c++)
        ClassAlive[c] = (c<nClusters2start);
}

// CEM(StartFile) - Does a whole CEM algorithm from a random start or masked start
// whereby clusters are assigned according to the similarity of their masks
// optional start file loads this cluster file to start iteration
// if Recurse is 0, it will not try and split.
// if InitRand is 0, use cluster assignments already in structure
scalar KK::CEM(char *CluFile, integer Recurse, integer InitRand,
        bool allow_assign_to_noise)
{
    integer p;
    integer nChanged;
    integer Iter;
    vector<integer> OldClass(nPoints);
    scalar Score, OldScore;
    integer LastStepFull; // stores whether the last step was a full one
    integer DidSplit;

    if (Debug)
    {
        Output("Entering CEM \n");
    }

    if (CluFile && *CluFile)
        LoadClu(CluFile);
    else if (InitRand)
    {
        // initialize data to random
        if((MaskStarts||UseMaskedInitialConditions) && (UseDistributional) && Recurse)
            StartingConditionsFromMasks();
        else
            StartingConditionsRandom();
    }

    // set all classes to alive
    Reindex();

    // main loop
    Iter = 0;
    FullStep = 1;
    Score = 0.0;
    do {
        // Store old classifications
        for(p=0; p<nPoints; p++) OldClass[p] = Class[p];

        // M-step - calculate class weights, means, and covariance matrices for each class
        MStep();

        // E-step - calculate scores for each point to belong to each class
        EStep();

        // dump distances if required

        if (DistDump) MatPrint(Distfp, &LogP.front(), DistDump, MaxPossibleClusters);

        // C-step - choose best class for each
        CStep(allow_assign_to_noise);

        // Compute class penalties
        ComputeClassPenalties();

        // Would deleting any classes improve things?
        if(Recurse) ConsiderDeletion();

        // Calculate number changed
        nChanged = 0;
        for(p=0; p<nPoints; p++) nChanged += (OldClass[p] != Class[p]);

        //Compute elapsed time
        timesofar = (clock()-Clock0)/(scalar) CLOCKS_PER_SEC;
        //Output("\nTime so far" SCALARFMT " seconds.\n", timesofar);
        
        //Write start of output to klg file
        if(Verbose>=1)
        {
            if(Recurse==0) Output("\t\tSP:");
            if ((Recurse!=0)||(SplitInfo==1&&Recurse==0))
                Output("Iteration %d%c (" SCALARFMT " sec): %d clusters\n",
				       (int)Iter, FullStep ? 'F' : 'Q', timesofar, (int)nClustersAlive);
        }
        
        // Calculate score
        OldScore = Score;
        Score = ComputeScore();
        
        //Finish output to klg file with Score already returned via the ComputeScore() function
        if(Verbose>=1)
        {
			Output("  nChanged %d\n", (int)nChanged);
        }

        //if(Verbose>=1)
        //{
        //    if(Recurse==0) Output("\t");
        //    Output(" Iteration %d%c: %d clusters Score %.7g nChanged %d\n",
        //        (int)Iter, FullStep ? 'F' : 'Q', (int)nClustersAlive, Score, (int)nChanged);
        //}

        Iter++;
        numiterations++;
        global_numiterations++;
        iteration_metric2 += (scalar)(nDims*nDims)*(scalar)(nPoints);
        iteration_metric3 += (scalar)(nDims*nDims)*(scalar)(nDims*nPoints);
        

        if (Debug)
        {
            for(p=0;p<nPoints;p++) BestClass[p] = Class[p];
            SaveOutput();
            Output("Press return");
            getchar();
        }

        // Next step a full step?
        LastStepFull = FullStep;
        FullStep = (
                        nChanged>ChangedThresh*nPoints
                        || nChanged == 0
                        || Iter%FullStepEvery==0
                        || Score > OldScore // SNK: Resurrected
                    //SNK    Score decreases ARE because of quick steps!
                    ) ;
        if (Iter>MaxIter)
        {
            Output("Maximum iterations exceeded\n");
            break;
        }
        
		//Save a temporary clu file when not splitting
		if (SaveTempCluEveryIter && Recurse)
		{
            SaveTempOutput(); //SNK Saves a temporary output clu file on each iteration
            Output("Writing temp and besttemp clu file \n");	
		}
		
        // try splitting
        //integer mod = (abs(Iter-SplitFirst))%SplitEvery;
        //Output("\n Iter mod SplitEvery = %d\n",(int)mod);
        //Output("Iter-SplitFirst %d \n",(int)(Iter-SplitFirst));
        if ((Recurse && SplitEvery>0) && ( Iter==SplitFirst  ||( Iter>=SplitFirst+1 && (Iter-SplitFirst)%SplitEvery==SplitEvery-1 )  || (nChanged==0 && LastStepFull) ) )
        {   
			SaveTempOutput(); //SNK Saves a temporary output clu file before each split
		    Output("Writing temp and besttemp clu file \n");	
            
            DidSplit = TrySplits();
        } else DidSplit = 0;

    } while (nChanged > 0 || !LastStepFull || DidSplit);

    if (DistDump) fprintf(Distfp, "\n");

    return Score;
}

// does the two-step clustering algorithm:
// first make a subset of the data, to SubPoints points
// then run CEM on this
// then use these clusters to do a CEM on the full data
// It calls CEM whenever there is no initialization clu file (i.e. the most common usage)
scalar KK::Cluster(char *StartCluFile=NULL)
{
    if (Debug)
    {
    Output("Entering Cluster \n");
    }
    
    
    
    if (Subset<=1)
    { // don't subset
		Output("------ Clustering full data set of %d points ------\n", (int)nPoints);
        return CEM(NULL, 1, 1);
    }

    // otherwise run on a subset of points
    integer sPoints = nPoints/Subset; // number of subset points - integer division will round down

    vector<integer> SubsetIndices(sPoints);
    for (integer i=0; i<sPoints; i++)
        // choose point to include, evenly spaced plus a random offset
        SubsetIndices[i] = Subset*i + irand(0, Subset-1);
    KK KKSub = KK(*this, SubsetIndices);

    // run CEM algorithm on KKSub
	Output("------ Running on subset of %d points ------\n", (int)sPoints);
    KKSub.CEM(NULL, 1, 1);

    // now copy cluster shapes from KKSub to main KK
    Weight = KKSub.Weight;
    Mean = KKSub.Mean;
    Cov = KKSub.Cov;
	DynamicCov = KKSub.DynamicCov;
    ClassAlive = KKSub.ClassAlive;
    nClustersAlive = KKSub.nClustersAlive;
    AliveIndex = KKSub.AliveIndex;

    // Run E and C steps on full data set
	Output("------ Evaluating fit on full set of %d points ------\n", (int)nPoints);
	if(UseDistributional)
		ComputeClusterMasks(); // needed by E-step normally computed by M-step
    EStep();
    CStep();

    // compute score on full data set and leave
    return ComputeScore();
}

// Initialise by loading data from files
KK::KK(char *FileBase, integer ElecNo, char *UseFeatures,
        scalar PenaltyK, scalar PenaltyKLogN, integer PriorPoint)
{
	KK_split = NULL;
	K2_container = NULL;
    penaltyK = PenaltyK;
    penaltyKLogN = PenaltyKLogN;
    LoadData(FileBase, ElecNo, UseFeatures);
    priorPoint = PriorPoint;
    
    //NOTE: penaltyK, penaltyKlogN, priorPoint, lower case versions of global variable PenaltyK PenaltyKLogN and PriorPoint
    
    DoInitialPrecomputations();//Now DoPrecomputations is only invoked in the initialization
    numiterations = 0;
	init_type = 0;
}

// This function is used by both of the constructors below, it initialises
// the data from a source KK object with a subset of the indices.
void KK::ConstructFrom(const KK &Source, const vector<integer> &Indices)
{
	KK_split = NULL;
	K2_container = NULL;
    nDims = Source.nDims;
    nDims2 = nDims*nDims;
    nPoints = Indices.size();
    penaltyK = Source.penaltyK;
    penaltyKLogN = Source.penaltyKLogN;
    priorPoint = Source.priorPoint;
    nStartingClusters = Source.nStartingClusters;
	NoisePoint = Source.NoisePoint;
	FullStep = Source.FullStep;
	nClustersAlive = Source.nClustersAlive;
	numiterations = Source.numiterations;
    AllocateArrays(); // Set storage for all the arrays such as Data, FloatMasks, Weight, Mean, Cov, etc.

    if (Debug)
    {
        Output("Entering ConstructFrom: \n");
    }

    // fill with a subset of points
    for (integer p=0; p<nPoints; p++)
    {
        integer psource = Indices[p];
        //copy data and masks
        for (integer d=0; d<nDims; d++)
            Data[p*nDims + d] = Source.Data[psource*nDims + d];
        for (integer d=0; d<nDims; d++)
            Masks[p*nDims + d] = Source.Masks[psource*nDims + d];
        if(UseDistributional)
        {
            for (integer d=0; d<nDims; d++)
            //    CorrectionTerm[p*nDims + d] = Source.CorrectionTerm[psource*nDims + d];
                FloatMasks[p*nDims + d] = Source.FloatMasks[psource*nDims + d];
        }
        
        UnMaskDims[p] = Source.UnMaskDims[psource];

    }

    

    //Output(" Printing Source.NoiseVariance[2] = %f",Source.NoiseVariance[2]);

    if(UseDistributional)
    {
        NoiseMean.resize(nDims);
        NoiseVariance.resize(nDims);
        nMasked.resize(nDims);

    
        for (integer d=0; d<nDims;d++)
        {
            NoiseMean[d] = Source.NoiseMean[d];
            NoiseVariance[d] = Source.NoiseVariance[d];
            nMasked[d] = Source.nMasked[d];

        }
    }

    DoPrecomputations();
    
    //Output(" Printing Source.NoiseMean[2] = %f",NoiseVariance[2]);

    numiterations = 0;
}

void KK::ConstructFrom(const KK &Source)
{
	vector<integer> Indices(Source.nPoints);
	for (integer i = 0; i<Source.nPoints; i++)
		Indices[i] = i;
	ConstructFrom(Source, Indices);
}

KK::KK(const KK &Source, const vector<integer> &Indices)
{
    ConstructFrom(Source, Indices);
	init_type = 2;
}

// If we don't specify an index subset, use everything.
KK::KK(const KK &Source)
{
	ConstructFrom(Source);
	init_type = 1;
}

KK::~KK()
{
	if (KK_split) delete KK_split;
	KK_split = NULL;
	if (K2_container) delete K2_container;
	K2_container = NULL;
}

// Main loop
int main(int argc, char **argv)
{
    scalar Score;
    scalar BestScore = HugeScore;
    integer p, i;
    SetupParams((integer)argc, argv); // This function is defined in parameters.cpp
	if (RamLimitGB == 0.0)
	{
		RamLimitGB = (1.0*available_physical_memory()) / (1024.0*1024.0*1024.0);
#ifdef __APPLE__
		Output("Setting RAM limit to total physical memory, %.2f GB.\n", (double)RamLimitGB);
		Output("WARNING: Not all physical memory will be available, but on Macs it is not possible\n");
		Output("         to get the available physical memory.\n");
#else
		Output("Setting RAM limit to available physical memory, %.2f GB.\n", (double)RamLimitGB);
#endif
	}
	else if (RamLimitGB < 0.0)
	{
		RamLimitGB = 1e20;
		Output("WARNING: You have chosen not to set a RAM limit, this may cause problems.\n");
	}
	memory_tracker.limit_gb = RamLimitGB;
    
    //clock_t Clock0 = clock();
    Clock0 = clock();
    

    // The main KK object, loads the data and does some precomputations
    KK K1(FileBase, ElecNo, UseFeatures, PenaltyK, PenaltyKLogN, PriorPoint);
    if(UseDistributional && SaveSorted) //Bug fix (Classical KK would terminate here)
        K1.SaveSortedData();

    // Seed random number generator
    srand((unsigned int)RandomSeed);

    // open distance dump file if required
    if (DistDump) Distfp = fopen("DISTDUMP", "w");

    // start with provided file, if required
    if (*StartCluFile)
    {
        Output("\nStarting from cluster file %s\n", StartCluFile);
        
        scalar iterationtime = (scalar)clock();
        BestScore = K1.CEM(StartCluFile, 1, 1);  //Main computation
        iterationtime = (clock()-iterationtime)/(scalar) CLOCKS_PER_SEC;
        Output("Time taken for this iteration:" SCALARFMT " seconds.\n", iterationtime);
        
		Output(" %d->%d Clusters: Score " SCALARFMT "\n\n", (int)K1.nStartingClusters, (int)K1.nClustersAlive, BestScore);
        for(p=0; p<K1.nPoints; p++)
            K1.BestClass[p] = K1.Class[p];
	    K1.SaveOutput();
    }
    else
	{
    // loop through numbers of clusters ...
    for(K1.nStartingClusters=(int)MinClusters; K1.nStartingClusters<=(int)MaxClusters; K1.nStartingClusters++)
        for(i=0; i<nStarts; i++)
        {
            // do CEM iteration
			Output("\nStarting from %d clusters...\n", (int)K1.nStartingClusters);
            scalar iterationtime = (scalar)clock();
            Score = K1.Cluster(); //Main computation
            iterationtime = (clock()-iterationtime)/(scalar) CLOCKS_PER_SEC;
            Output("Time taken for this iteration:" SCALARFMT " seconds.\n", iterationtime);

			Output(" %d->%d Clusters: Score " SCALARFMT ", best is " SCALARFMT "\n", (int)K1.nStartingClusters, (int)K1.nClustersAlive, Score, BestScore);
            if (Score < BestScore)
            {
                Output("THE BEST YET!\n"); // New best classification found
                BestScore = Score;
                for(p=0; p<K1.nPoints; p++)
                    K1.BestClass[p] = K1.Class[p];
                K1.SaveOutput();
            }
            Output("\n");
        }
	}	

    K1.SaveOutput();

    scalar tottime = (clock()-Clock0)/(scalar) CLOCKS_PER_SEC;

    Output("Main iterations: %d (time per iteration =" SCALARFMT " ms)\n",
		    (int)K1.numiterations,
            1e3*tottime/K1.numiterations);
    Output("Total iterations: %d (time per iteration =" SCALARFMT " ms)\n",
		    (int)global_numiterations,
            1e3*tottime/global_numiterations);
    Output("\nDef. Iteration metric 2:\nIteration_metric2 += (scalar)(nDims*nDims)*(scalar)(nPoints)\n");
    Output("Iterations metric 2: " SCALARFMT " (time per metric unit =" SCALARFMT "ns)\n",
            iteration_metric2,
            1e9*tottime/iteration_metric2);
     Output("\nDef. Iteration metric 3:\nIteration_metric3 += (scalar)(nDims*nDims)*(scalar)(nDims*nPoints)\n");
    Output("Iterations metric 3: " SCALARFMT " (time per metric unit=" SCALARFMT "ps)\n",
            iteration_metric3,
            1e12*tottime/iteration_metric3);
    Output("\nThat took " SCALARFMT " seconds.\n", tottime);

    if (DistDump) fclose(Distfp);

    return 0;
}
