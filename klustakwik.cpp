// MaskedKlustaKwik2.C
//
// Fast clustering using the CEM algorithm with Masks.

#include "klustakwik.h"

// GLOBAL VARIABLES
FILE *Distfp;
int global_numiterations = 0;
scalar iteration_metric2 = (scalar)0;
scalar iteration_metric3 = (scalar)0;

// Sets storage for KK class.  Needs to have nDims and nPoints defined
void KK::AllocateArrays() {

	nDims2 = nDims*nDims;
	NoisePoint = 1; // Ensures that the mixture weight for the noise cluster never gets to zero

	// Set sizes for arrays
	Data.resize(nPoints * nDims);
	//SNK
	Masks.resize(nPoints * nDims);
	FloatMasks.resize(nPoints * nDims);
	UnMaskDims.resize(nPoints); //SNK Number of unmasked dimensions for each data point when using float masks $\sum m_i$
	Weight.resize(MaxPossibleClusters);
	Mean.resize(MaxPossibleClusters*nDims);
	Cov.resize(MaxPossibleClusters*nDims2);
	LogP.resize(MaxPossibleClusters*nPoints);
	Class.resize(nPoints);
	OldClass.resize(nPoints);
	Class2.resize(nPoints);
	BestClass.resize(nPoints);
	ClassAlive.resize(MaxPossibleClusters);
	AliveIndex.resize(MaxPossibleClusters);
	ClassPenalty.resize(MaxPossibleClusters);
	nClassMembers.resize(MaxPossibleClusters);
	if(UseDistributional)
	{
		CorrectionTerm.resize(nPoints * nDims);
		ClassCorrectionFactor.resize(MaxPossibleClusters*nDims);
	}
}

// recompute index of alive clusters (including 0, the noise cluster)
// should be called after anything that changes ClassAlive
void KK::Reindex()
{
    int c;

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
scalar KK::Penalty(int n)
{
	int nParams;
	if(n==1)
		return 0;
	nParams = (nDims*(nDims+1)/2 + nDims + 1)*(n-1); // each has cov, mean, &p 

	scalar p = penaltyK*(scalar)(nParams*2) // AIC
			+penaltyKLogN*((scalar)nParams*(scalar)log((scalar)nPoints)/2); // BIC
	return p;
}


// Penalties for Masked CEM
void KK::ComputeClassPenalties()
{
	if(!((bool)UseDistributional)) // This function must only be called in Use Distributional  mode
		return;
	for(int c=0; c<MaxPossibleClusters; c++)
		ClassPenalty[c] = (scalar)0;
	// compute sum of nParams for each
	vector<int> NumberInClass(MaxPossibleClusters);
	for(int p=0; p<nPoints; p++)
	{
		int c = Class[p];
		NumberInClass[c]++;
	//	int n = UnmaskedInd[p+1]-UnmaskedInd[p]; // num unmasked dimensions
        scalar n = UnMaskDims[p];
		scalar nParams = n*(n+1)/2+n+1;
		ClassPenalty[c] += nParams;
	}
	// compute mean nParams for each cluster
	for(int c=0; c<MaxPossibleClusters; c++)
		if(NumberInClass[c]>0)
			ClassPenalty[c] /= (scalar)NumberInClass[c];
	// compute penalty for each cluster
	for(int c=0; c<MaxPossibleClusters; c++)
	{
		scalar nParams = ClassPenalty[c];
		ClassPenalty[c] = penaltyK*(scalar)(nParams*2)
				+penaltyKLogN*((scalar)nParams*(scalar)log((scalar)nPoints)/2);
	}
}




// M-step: Calculate mean, cov, and weight for each living class
// also deletes any classes with fewer points than nDim
void KK::MStep()
{
	int p, c, cc, i, j;
	vector<scalar> Vec2Mean(nDims);

	// clear arrays
	memset((void*)&nClassMembers.front(), 0, MaxPossibleClusters*sizeof(int));
	memset((void*)&Mean.front(), 0, MaxPossibleClusters*nDims*sizeof(scalar));
	memset((void*)&Cov.front(), 0, MaxPossibleClusters*nDims*nDims*sizeof(scalar));
// NOTE: memset commands above replace the code below:
//	for(c=0; c<MaxPossibleClusters; c++) {
//		nClassMembers[c] = 0;
//		for(i=0; i<nDims; i++) Mean[c*nDims + i] = 0;
//		for(i=0; i<nDims; i++) for(j=i; j<nDims; j++) {
//			Cov[c*nDims2 + i*nDims + j] = 0;
//		}
//	}

	if (Debug) { Output("Entering Unmasked Mstep \n");}

	// Accumulate total number of points in each class
	for (p=0; p<nPoints; p++) nClassMembers[Class[p]]++;

    // check for any dead classes

    if(UseDistributional)
    {
    	for (int cc=0; cc<nClustersAlive; cc++)
        {
            int c = AliveIndex[cc];
            if (Debug){Output("DistributionalMstep: Class %d contains %d members \n", c, nClassMembers[c]);}
            	if (c>0 && nClassMembers[c]<1)//nDims)
            	{
            		ClassAlive[c]=0;
            		if (Debug) {Output("UnmaskedMstep_dist: Deleted class %d: no members\n", c);}
            	}
        }
    }
    else
    {

    	for (cc=0; cc<nClustersAlive; cc++)
        {
    	        c = AliveIndex[cc];
    	        if (Debug) {Output("Mstep: Class %d contains %d members \n", c, nClassMembers[c]);}
    	        if (c>0 && nClassMembers[c]<=nDims)
                {
    	            ClassAlive[c]=0;
    	        	if (Debug) {Output("Deleted class %d: not enough members\n", c);}
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
//	for (p=0; p<nPoints; p++)
//	{
//		c = Class[p];
//		// calculate distance from mean
//		for(i=0; i<nDims; i++)
//			Vec2Mean[i] = Data[p*nDims + i] - Mean[c*nDims + i];
//		for(i=0; i<nDims; i++)
//			for(j=i; j<nDims; j++)
//				Cov[c*nDims2 + i*nDims + j] += Vec2Mean[i] * Vec2Mean[j];
//	}
    if((int)AllVector2Mean.size()<nPoints*nDims)
    	AllVector2Mean.resize(nPoints*nDims);
    vector< vector<int> > PointsInClass(MaxPossibleClusters);
    for(p=0; p<nPoints; p++)
    {
    	c = Class[p];
    	PointsInClass[c].push_back(p);
		for(i=0; i<nDims; i++)
			AllVector2Mean[p*nDims+i] = Data[p*nDims + i] - Mean[c*nDims + i];
    }
    for(c=0; c<MaxPossibleClusters; c++)
    {
    	vector<int> &PointsInThisClass = PointsInClass[c];
    	SafeArray<scalar> safeCov(Cov, c*nDims2, "safeCovMStep");
    	for(int iblock=0; iblock<nDims; iblock+=COVARIANCE_BLOCKSIZE)
    		for(int jblock=iblock; jblock<nDims; jblock+=COVARIANCE_BLOCKSIZE)
				for(int q=0; q<(int)PointsInThisClass.size(); q++)
				{
					p = PointsInThisClass[q];
					scalar *cv2m = &AllVector2Mean[p*nDims];
					for(i=iblock; i<MIN(nDims, iblock+COVARIANCE_BLOCKSIZE); i++)
					{
						scalar cv2mi = cv2m[i];
						int jstart;
						if(jblock!=iblock)
							jstart = jblock;
						else
							jstart = i;
						scalar *covptr = &safeCov[i*nDims+jstart];
						scalar *cv2mjptr = &cv2m[jstart];
						//scalar *cv2mjend = cv2m+MIN(nDims, jblock+COVARIANCE_BLOCKSIZE);
						//for(j=jstart; j<MIN(nDims, jblock+COVARIANCE_BLOCKSIZE); j++)
						//for(; cv2mjptr!=cv2mjend;)
						for(j=MIN(nDims, jblock+COVARIANCE_BLOCKSIZE)-jstart; j; j--)
							*covptr++ += cv2mi*(*cv2mjptr++);
					}
				}
    }

    if(UseDistributional)
    {
    	for(cc=0; cc<nClustersAlive; cc++)
    	{
    		c = AliveIndex[cc];
        	vector<int> &PointsInThisClass = PointsInClass[c];
    		for(i=0; i<nDims; i++)
    		{
    			scalar ccf = 0.0; // class correction factor
				for(int q=0; q<(int)PointsInThisClass.size(); q++)
				{
					p = PointsInThisClass[q];
					ccf += CorrectionTerm[p*nDims+i];
				}
				//Output("Class %d Class correction factor[%d] = %f \n",c,i,ccf);
				Cov[c*nDims2+i*nDims+i] += ccf;
			//	Output("Class %d Covariance diagonal[%d] = %f \n",c,i,Cov[c*nDims2+i*nDims+i] );
                ClassCorrectionFactor[c*nDims+i] = ccf/(scalar)(nClassMembers[c]*nClassMembers[c]);
    		}
    	}

    // Add a diagonal matrix of Noise variances to the covariance matrix for renormalization
    	for (cc=0; cc<nClustersAlive; cc++)
    	        {c = AliveIndex[cc];
    	    	for (i=0; i<nDims; i++)
    	    		{
    	    		//Output("Class %d: PriorPoint*NoiseVariance[%d] = %f",c,i,priorPoint*NoiseVariance[i]);
    	    		Cov[c*nDims2+i*nDims+i] += priorPoint*NoiseVariance[i];
    	    		}
    	        }


    }

    // and normalize
    if(UseDistributional){
        	    for (cc=0; cc<nClustersAlive; cc++)
        	    {
        	        c = AliveIndex[cc];
        			for(i=0; i<nDims; i++)
        				for(j=i; j<nDims; j++)
        					Cov[c*nDims2 + i*nDims + j] /= (nClassMembers[c]+priorPoint-1);
        		}

    	}
    else {    //For original KlustaKwik classical EM
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
			Output("Class %d - Weight %.2g\n", c, Weight[c]);
			Output("Mean: ");
			MatPrint(stdout, &Mean.front() + c*nDims, 1, nDims);
			Output("\nCov:\n");
			MatPrint(stdout, &Cov.front() + c*nDims2, nDims, nDims);
			Output("\n");
		}
	}
}



// E-step.  Calculate Log Probs for each point to belong to each living class
// will delete a class if covariance matrix is singular
// also counts number of living classes
void KK::EStep()
{
	int p, c, cc, i;
	int nSkipped;
	scalar LogRootDet; // log of square root of covariance determinant
	scalar Mahal; // Mahalanobis distance of point from cluster center
	scalar correction_factor = (scalar)1; // for partial correction in distributional step
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
	for (p=0; p<nPoints; p++)
		LogP[p*MaxPossibleClusters + 0] = (float)-log(Weight[0]);

    for(cc=1; cc<nClustersAlive; cc++)
    {
        c = AliveIndex[cc];

		// calculate cholesky decomposition for class c
        SafeArray<scalar> safeCov(Cov, c*nDims2, "safeCov");
		if(Cholesky(safeCov, safeChol, nDims))
		{
			// If Cholesky returns 1, it means the matrix is not positive definite.
			// So kill the class.
            // Cholesky is defined in linalg.cpp
			Output("Unmasked E-step: Deleting class %d: covariance matrix is singular \n", c);
			ClassAlive[c] = 0;
			continue;
		}

		// LogRootDet is given by log of product of diagonal elements
		LogRootDet = 0;
		for(i=0; i<nDims; i++)
			LogRootDet += (float)log(Chol[i*nDims + i]);

		// if distributional E step, compute diagonal of inverse of cov matrix
		if(UseDistributional)
		{
			vector<scalar> BasisVector(nDims);
			SafeArray<scalar> safeBasisVector(BasisVector, "BasisVector");
			for(int i=0; i<nDims; i++)
				safeBasisVector[i] = (scalar)0;
			for(int i=0; i<nDims; i++)
			{   safeBasisVector[i] = (scalar)1;
				// calculate Root vector - by Chol*Root = BasisVector
				TriSolve(safeChol, safeBasisVector, safeRoot, nDims);
				// add half of Root vector squared to log p
				scalar Sii = (scalar)0;
				for(int j=0; j<nDims; j++)
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

			// Compute Mahalanobis distance
			Mahal = 0;

			// calculate data minus class mean
			for(i=0; i<nDims; i++)
				Vec2Mean[i] = Data[p*nDims + i] - Mean[c*nDims + i];

			// calculate Root vector - by Chol*Root = Vec2Mean
			TriSolve(safeChol, safeVec2Mean, safeRoot, nDims);

			// add half of Root vector squared to log p
			for(i=0; i<nDims; i++)
				Mahal += Root[i]*Root[i];
    //        if(Debug)Output("Mahal = %f",Mahal);

			// if distributional E step, add correction term
			if(UseDistributional)
				for(i=0; i<nDims; i++)
				{
					//if(UseDistributionalEStep==2)      // Distribution E-Step 2 "semi-Bayesian", no longer used, code retained here in case
					//	correction_factor = ClassCorrectionFactor[c*nDims+i]+
					//			(1.0-2.0/(scalar)nClassMembers[c]);
					Mahal += correction_factor*CorrectionTerm[p*nDims+i]*safeInvCovDiag[i];
		//			                    if(Debug) {Output("CorrectionTerm[%d*nDims+%d] = %f ",p,CorrectionTerm[p*nDims+i],i);
		//			   Output("Mahal = %f",Mahal);}
                }

			// Score is given by Mahal/2 + log RootDet - log weight
			LogP[p*MaxPossibleClusters + c] = Mahal/2
   									+ LogRootDet
									- log(Weight[c])
									+ (float)log(2*M_PI)*nDims/2;
			                              //           Output("LogP = %d ",LogP[p*MaxPossibleClusters + c]);

		} // for(p=0; p<nPoints; p++)
	} // for(cc=1; cc<nClustersAlive; cc++)
}



// Choose best class for each point (and second best) out of those living
void KK::CStep(bool allow_assign_to_noise)
{
	int p, c, cc, TopClass, SecondClass;
	int ccstart = 0;
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

	int c, p, CandidateClass=0;
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
	for(p=0; p<nPoints; p++)
    {
		DeletionLoss[Class[p]] += LogP[p*MaxPossibleClusters + Class2[p]] - LogP[p*MaxPossibleClusters + Class[p]];
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

	//Output("cand Class %d would lose " SCALARFMT " gain is " SCALARFMT "\n", CandidateClass, Loss, DeltaPen);
	// is it worth it?
    //06/12/12 fixing bug introduced which considered DeltaPen twice!
	if (UseDistributional) //For the distributional algorithm we need to use the ClusterPenalty
    {
        if (Loss<0)
        {
            Output("Deleting Class %d. Lose " SCALARFMT " but Gain " SCALARFMT "\n", CandidateClass, DeletionLoss[CandidateClass], DeltaPen);
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
            Output("Deleting Class %d. Lose " SCALARFMT " but Gain " SCALARFMT "\n", CandidateClass, DeletionLoss[CandidateClass], DeltaPen);
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
    int p, c, val;
    int status;


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
int KK::TrySplits()
{
    int c, cc, c2, p, p2, DidSplit = 0;
    scalar Score, NewScore, UnsplitScore, SplitScore;
    int UnusedCluster;
    //KK K2; // second KK structure for sub-clustering
    //KK K3; // third one for comparison

    if(nClustersAlive>=MaxPossibleClusters-1)
    {
        Output("Won't try splitting - already at maximum number of clusters\n");
        return 0;
    }

    // set up K3 and remember to add the masks
    KK K3(*this);
    Output("Compute initial score before splitting: ");
    Score = ComputeScore();

    // loop thu clusters, trying to split
    for (cc=1; cc<nClustersAlive; cc++)
    {
        c = AliveIndex[cc];

        // set up K2 structure to contain points of this cluster only

        vector<int> SubsetIndices;
        for(p=0; p<nPoints; p++)
        	if(Class[p]==c)
        		SubsetIndices.push_back(p);
        if(SubsetIndices.size()==0)
        	continue;
        KK K2(*this, SubsetIndices);

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
        if (Verbose>=1) Output("\n Trying to split cluster %d (%d points) \n", c, K2.nPoints);
        K2.nStartingClusters=2; // (2 = 1 clusters + 1 unused noise cluster)
        UnsplitScore = K2.CEM(NULL, 0, 1, false);
        K2.nStartingClusters=3; // (3 = 2 clusters + 1 unused noise cluster)
        SplitScore = K2.CEM(NULL, 0, 1, false);

        // Fix by Michaël Zugaro: replace next line with following two lines
        // if(SplitScore<UnsplitScore) {
        if(K2.nClustersAlive<2) Output("\n Split failed - leaving alone\n");
        if((SplitScore<UnsplitScore)&&(K2.nClustersAlive>=2)) {
            // will splitting improve the score in the whole data set?

            // assign clusters to K3
            for(c2=0; c2<MaxPossibleClusters; c2++) K3.ClassAlive[c2]=0;
         //   Output("%d Points in class %d in KKobject K3 ", c2, K3.nClassMembers[c2]);
            p2 = 0;
            for(p=0; p<nPoints; p++)
            {
                if(Class[p]==c)
                {
                    if(K2.Class[p2]==1) K3.Class[p] = c;
                    else if(K2.Class[p2]==2) K3.Class[p] = UnusedCluster;
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
            Output("About to compute K3 class penalties");
            if (UseDistributional) K3.ComputeClassPenalties(); //SNK Fixed bug: Need to compute the cluster penalty properly, cluster penalty is only used in UseDistributional mode
            NewScore = K3.ComputeScore();
            Output("\n Splitting cluster %d changes total score from " SCALARFMT " to " SCALARFMT "\n", c, Score, NewScore);

            if (NewScore<Score)
            {
                DidSplit = 1;
                Output("\n So it's getting split into cluster %d.\n", UnusedCluster);
                // so put clusters from K3 back into main KK struct (K1)
                for(c2=0; c2<MaxPossibleClusters; c2++) ClassAlive[c2] = K3.ClassAlive[c2];
                for(p=0; p<nPoints; p++) Class[p] = K3.Class[p];
            } else
            {
                Output("\n So it's not getting split.\n");
            }
        }
    }
    return DidSplit;
}

// ComputeScore() - computes total score.  Requires M, E, and C steps to have been run
scalar KK::ComputeScore()
{
    int p;
   // int debugadd;

    scalar penalty = (scalar)0;
    if(UseDistributional)  // For distributional algorithm we require the cluster penalty
		for(int c=0; c<MaxPossibleClusters; c++)
			penalty += ClassPenalty[c];
    else
    	penalty = Penalty(nClustersAlive);
    scalar Score = penalty;
    for(p=0; p<nPoints; p++)
    {	//debugadd = LogP[p*MaxPossibleClusters + Class[p]];
        Score += LogP[p*MaxPossibleClusters + Class[p]];
		// Output("point %d: cumulative score " SCALARFMT " adding" SCALARFMT "\n", p, Score, debugadd);
    }
    //Error("Score: " SCALARFMT " Penalty: " SCALARFMT "\n", Score, penalty);
    Output("Score: Raw " SCALARFMT " + Penalty " SCALARFMT " = " SCALARFMT, Score-penalty, penalty, Score);

	if (Debug) {
		int c, cc;
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
	    for(int p=0; p<nPoints; p++) // No points are put in the noise cluster to begin 
	    	Class[p] = irand(1, nStartingClusters-1);
    else
        for(int p=0; p<nPoints; p++) //If there is only one cluster, put all the points in the noise cluster
        	Class[p] = 0;

	for(int c=0; c<MaxPossibleClusters; c++)
		ClassAlive[c] = (c<nStartingClusters);

	Output("Assigned %d initial classes randomly.\n", nStartingClusters);
}

// Initialise starting conditions by selecting unique masks at random
void KK::StartingConditionsFromMasks()
{
    int nClusters2start; //SNK To replace nStartingClusters within this variable only
    
    //if (Debug)
    //    Output("StartingConditionsFromMasks: ");
    Output("Starting initial clusters from distinct float masks \n ");
    
	if(nStartingClusters<=1) // If only 1 starting clutser has been requested, assign all the points to cluster 0
	{
        for(int p=0; p<nPoints; p++)
        	Class[p] = 0;
	}
	else   
	{
		int num_masks = 0;
		for(int p=0; p<nPoints; p++)
			num_masks += (int)SortedMaskChange[p];
		
        if((nStartingClusters-1)>num_masks)
		{
			Error("Not enough masks (%d) to generate starting clusters (%d), "
					"so starting with (%d) clusters instead.\n", num_masks,
					nStartingClusters, num_masks+1);
            nClusters2start = num_masks+1;
			//return;
		}
        else
        {
            nClusters2start = nStartingClusters;
        }
		// Construct the set of all masks
		vector<bool> MaskUsed;
		vector<int> MaskIndex(nPoints);
		vector<int> MaskPointIndex;
		int current_mask_index = -1;
		for(int q=0; q<nPoints; q++)
		{
			int p = SortedIndices[q];
			if(q==0 || SortedMaskChange[p])
			{
				current_mask_index++;
				MaskUsed.push_back(false);
				MaskPointIndex.push_back(p);
			}
			MaskIndex[p] = current_mask_index;
		}
		// Select points at random until we have enough masks
		int masks_found = 0;
		vector<int> MaskIndexToUse;
		vector<int> FoundMaskIndex(num_masks);
		while(masks_found<nClusters2start-1)
		{
			int p = irand(0, nPoints-1);
			int mask_index = MaskIndex[p];
			if(!MaskUsed[mask_index])
			{
				MaskIndexToUse.push_back(mask_index);
				MaskUsed[mask_index] = true;
				FoundMaskIndex[mask_index] = masks_found;
				masks_found++;
			}
		}
		// Assign points to clusters based on masks
		for(int p=0; p<nPoints; p++)
		{
			if(MaskUsed[MaskIndex[p]]) // we included this points mask
				Class[p] = FoundMaskIndex[MaskIndex[p]]+1; // so assign class to mask index
			else // this points mask not included
			{
				// so find closest match
				int closest_index = 0;
				int distance = nDims+1;
				vector<int> possibilities;
				for(int mi=0; mi<nClusters2start-1; mi++)
				{
					int mip = MaskPointIndex[MaskIndexToUse[mi]];
					// compute mask distance
					int curdistance = 0;
					for(int i=0; i<nDims; i++)
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
				if(AssignToFirstClosestMask)
					closest_index = possibilities[0];
				else
					closest_index = possibilities[irand(0, possibilities.size()-1)];
				Class[p] = closest_index+1;
			}
		}
		// print some info
		Output("Assigned %d initial classes from %d unique masks.\n",
				nClusters2start, num_masks);
		// Dump initial random classes to a file - knowledge of maskstart configuration may be useful
		// TODO: remove this for final version - SNK: actually it is a nice idea to keep this
		char fname[STRLEN];
		FILE *fp;
		sprintf(fname, "%s.initialclusters.%d.clu.%d", FileBase, nClusters2start, ElecNo);
		fp = fopen_safe(fname, "w");
		fprintf(fp, "%d\n", nClusters2start);
		for(int p=0; p<nPoints; p++)
			fprintf(fp, "%d\n", Class[p]);
		fclose(fp);
	}
	for(int c=0; c<MaxPossibleClusters; c++)
		ClassAlive[c] = (c<nClusters2start);
}

// CEM(StartFile) - Does a whole CEM algorithm from a random start or masked start
// whereby clusters are assigned according to the similarity of their masks
// optional start file loads this cluster file to start iteration
// if Recurse is 0, it will not try and split.
// if InitRand is 0, use cluster assignments already in structure
scalar KK::CEM(char *CluFile, int Recurse, int InitRand,
		bool allow_assign_to_noise)
{
	int p;
	int nChanged;
	int Iter;
	vector<int> OldClass(nPoints);
	scalar Score, OldScore;
	int LastStepFull; // stores whether the last step was a full one
    int DidSplit;

    if (Debug)
    {
        Output("Entering CEM \n");
    }

    if (CluFile && *CluFile)
    	LoadClu(CluFile);
	else if (InitRand)
	{
        // initialize data to random
		if(UseMaskedInitialConditions && (UseDistributional) && Recurse)
		//if(UseMaskedInitialConditions && (UseMaskedMStep||UseMaskedEStep) && Recurse)
			StartingConditionsFromMasks();
		else
			StartingConditionsRandom();
    }

	// set all classes to alive
    Reindex();

	// main loop
	Iter = 0;
	FullStep = 1;
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

        
        //Write start of output to klg file
        if(Verbose>=1)
        {
            if(Recurse==0) Output("\t");
            Output(" Iteration %d%c: %d clusters ",
                   Iter, FullStep ? 'F' : 'Q', nClustersAlive);
        }
        
		// Calculate score
		OldScore = Score;
		Score = ComputeScore();
        
        //Finish output to klg file with Score already returned via the ComputeScore() function
        if(Verbose>=1)
        {
            Output(" nChanged %d\n", nChanged);
        }

		//if(Verbose>=1)
        //{
        //    if(Recurse==0) Output("\t");
        //    Output(" Iteration %d%c: %d clusters Score %.7g nChanged %d\n",
		//	    Iter, FullStep ? 'F' : 'Q', nClustersAlive, Score, nChanged);
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
					//SNK	Score decreases ARE because of quick steps!
					) ;
		if (Iter>MaxIter)
        {
			Output("Maximum iterations exceeded\n");
			break;
		}

        // try splitting
        if ((Recurse && SplitEvery>0) && (Iter%SplitEvery==SplitEvery-1 || (nChanged==0 && LastStepFull)))
        {
            SaveTempOutput(); //SNK Saves a temporary output clu file before each split
            Output("Writing temp clu file \n");
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
		Output("--- Clustering full data set of %d points ---\n", nPoints);
		return CEM(NULL, 1, 1);
	}

	// otherwise run on a subset of points
	int sPoints = nPoints/Subset; // number of subset points - integer division will round down

	vector<int> SubsetIndices(sPoints);
	for (int i=0; i<sPoints; i++)
		// choose point to include, evenly spaced plus a random offset
		SubsetIndices[i] = Subset*i + irand(0, Subset-1);
	KK KKSub = KK(*this, SubsetIndices);

	// run CEM algorithm on KKSub
	Output("--- Running on subset of %d points ---\n", sPoints);
	KKSub.CEM(NULL, 1, 1);

	// now copy cluster shapes from KKSub to main KK
	Weight = KKSub.Weight;
	Mean = KKSub.Mean;
	Cov = KKSub.Cov;
	ClassAlive = KKSub.ClassAlive;
	nClustersAlive = KKSub.nClustersAlive;
	AliveIndex = KKSub.AliveIndex;

	// Run E and C steps on full data set
	Output("--- Evaluating fit on full set of %d points ---\n", nPoints);
	EStep();
	CStep();

	// compute score on full data set and leave
	return ComputeScore();
}

// Initialise by loading data from files
KK::KK(char *FileBase, int ElecNo, char *UseFeatures,
		scalar PenaltyK, scalar PenaltyKLogN, int PriorPoint)
{
	penaltyK = PenaltyK;
	penaltyKLogN = PenaltyKLogN;
	LoadData(FileBase, ElecNo, UseFeatures);
	priorPoint = PriorPoint;
	
	DoInitialPrecomputations();//Now DoPrecomputations is only invoked in the initialization
	numiterations = 0;
}

// This function is used by both of the constructors below, it initialises
// the data from a source KK object with a subset of the indices.
void KK::ConstructFrom(const KK &Source, const vector<int> &Indices)
{
	
	nDims = Source.nDims;
	nPoints = Indices.size();
	penaltyK = Source.penaltyK;
	penaltyKLogN = Source.penaltyKLogN;
	priorPoint = Source.priorPoint;
	nStartingClusters = Source.nStartingClusters;
	AllocateArrays(); // Set storage for all the arrays such as Data, FloatMasks, Weight, Mean, Cov, etc.

	if (Debug)
    {
        Output("Entering ConstructFrom: \n");
    }

	// fill with a subset of points
	for (int p=0; p<nPoints; p++)
	{
		int psource = Indices[p];
		//copy data and masks
		for (int d=0; d<nDims; d++)
			Data[p*nDims + d] = Source.Data[psource*nDims + d];
		for (int d=0; d<nDims; d++)
			Masks[p*nDims + d] = Source.Masks[psource*nDims + d];
		if(UseDistributional)
		{
			for (int d=0; d<nDims; d++)
			//	CorrectionTerm[p*nDims + d] = Source.CorrectionTerm[psource*nDims + d];
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

	
		for (int d=0; d<nDims;d++)
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

KK::KK(const KK &Source, const vector<int> &Indices)
{
	ConstructFrom(Source, Indices);
}

// If we don't specify an index subset, use everything.
KK::KK(const KK &Source)
{
	vector<int> Indices(Source.nPoints);
	for(int i=0; i<Source.nPoints; i++)
		Indices[i] = i;
	ConstructFrom(Source, Indices);
}

// Main loop
int main(int argc, char **argv)
{
	scalar Score;
	scalar BestScore = HugeScore;
	int p, i;
	SetupParams(argc, argv); // This function is defined in parameters.cpp
	
	clock_t Clock0 = clock();

	// The main KK object, loads the data and does some precomputations
	KK K1(FileBase, ElecNo, UseFeatures, PenaltyK, PenaltyKLogN, PriorPoint);
	if(SaveSorted)
		K1.SaveSortedData();

	// Seed random number generator
	srand(RandomSeed);

	// open distance dump file if required
	if (DistDump) Distfp = fopen("DISTDUMP", "w");

    // start with provided file, if required
    if (*StartCluFile)
    {
        Output("Starting from cluster file %s\n", StartCluFile);
        BestScore = K1.CEM(StartCluFile, 1, 1);
		Output(" %d->%d Clusters: Score " SCALARFMT "\n\n", K1.nStartingClusters, K1.nClustersAlive, BestScore);
		for(p=0; p<K1.nPoints; p++)
            K1.BestClass[p] = K1.Class[p];
		K1.SaveOutput();
    }

	// loop through numbers of clusters ...
	for(K1.nStartingClusters=MinClusters; K1.nStartingClusters<=MaxClusters; K1.nStartingClusters++)
        for(i=0; i<nStarts; i++)
        {
            // do CEM iteration
            Output("Starting from %d clusters...\n", K1.nStartingClusters);
		    Score = K1.Cluster();

		    Output(" %d->%d Clusters: Score " SCALARFMT ", best is " SCALARFMT "\n", K1.nStartingClusters, K1.nClustersAlive, Score, BestScore);
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

	K1.SaveOutput();

	scalar tottime = (clock()-Clock0)/(scalar) CLOCKS_PER_SEC;

	Output("Main iterations: %d (time per=" SCALARFMT "ms)\n",
			K1.numiterations,
			1e3*tottime/K1.numiterations);
	Output("Total iterations: %d (time per=" SCALARFMT "ms)\n",
			global_numiterations,
			1e3*tottime/global_numiterations);
	Output("Iterations metric 2: " SCALARFMT " (time per=" SCALARFMT "ns)\n",
			iteration_metric2,
			1e9*tottime/iteration_metric2);
	Output("Iterations metric 3: " SCALARFMT " (time per=" SCALARFMT "ps)\n",
			iteration_metric3,
			1e12*tottime/iteration_metric3);
	Output("\nThat took " SCALARFMT " seconds.\n", tottime);

	if (DistDump) fclose(Distfp);

	return 0;
}
