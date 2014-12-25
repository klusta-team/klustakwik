/*
 * io.cpp
 *
 * Handles input and output to files.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

// Disable some Visual Studio warnings
#define _CRT_SECURE_NO_WARNINGS

#include "klustakwik.h"
#include "numerics.h"

unsigned char convert_to_char(scalar x)
{
	integer y = (integer)(x*255.0);
	if(y<0) y = 0;
	if(y>255) y= 255;
	return (unsigned char)y;
}

// Loads in Fet file.  Also allocates storage for other arrays
void KK::LoadData(char *FileBase, integer ElecNo, char *UseFeatures)
{
    char fname[STRLEN];
    char fnamefmask[STRLEN];
    char line[STRLEN];
    integer p, i, j;
	// nFeatures is read as a %d so it has to be int type, not integer type
    int nFeatures, nmaskFeatures; // not the same as nDims! we don't use all features.
    FILE *fp;
    FILE *fpfmask;
    integer status;
    scalar val;
    integer UseLen;

    // open file
    sprintf(fname, "%s.fet.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "r");
    
    if((MaskStarts > 0)&& UseDistributional)
    {
        Output("-------------------------------------------------------------------------");
        Output("\nUsing Distributional EM with Maskstarts\n");
        MinClusters = MaskStarts;
        MaxClusters = MaskStarts;
        Output("NOTE: Maskstarts overides above values of MinClusters and MaxClusters \
                \nMinClusters = %d \nMaxClusters = %d \n ", (int)MinClusters, (int)MaxClusters);
    }
    
    if(UseDistributional)// replaces if(UseFloatMasks)
    {
        sprintf(fnamefmask,"%s.fmask.%d", FileBase, (int)ElecNo);
        fpfmask = fopen_safe(fnamefmask, "r");
    }
    else
    {
        fpfmask = NULL;
    }

    // count lines;
    nPoints=-1; // subtract 1 because first line is number of features
    while(fgets(line, STRLEN, fp)) {
        nPoints++;
    }

    // rewind file
    fseek(fp, 0, SEEK_SET);

    // read in number of features
    fscanf(fp, "%d", &nFeatures);
    if(Debug) Output("Number of features read in: %d \n ", nFeatures);

    // calculate number of dimensions
    if (UseFeatures[0] == 0)
    {
        nDims = nFeatures-DropLastNFeatures; // Use all but the last N Features.
        UseLen = nFeatures-DropLastNFeatures;
       // Output("nDims = %d ,UseLen = %d ", (int)nDims, (int)UseLen);
      //  UseFeatures =
    }
    else
    {
        UseLen = strlen(UseFeatures);
        nDims=0;
        for(i=0; i<nFeatures; i++)
        {
            nDims += (i<UseLen && UseFeatures[i]=='1');
        }
      //  Output("nDims = %d ,UseLen = %d ", (int)nDims, (int)UseLen);
    }
    nDims2 = nDims*nDims;
	MemoryCheck();
    AllocateArrays();

#ifdef STORE_DATA_AS_INTEGER
	// we need to scan through the data to find the min and max of each dimension before we save to memory
	vector<scalar> dmin(nFeatures);
	vector<scalar> dmax(nFeatures);
	for(p=0; p<nPoints; p++)
	{
		for(i=0; i<nFeatures; i++)
		{
            float readfloatval;
            status = fscanf(fp, "%f", &readfloatval);
            val = (scalar)readfloatval;
            if (status==EOF) Error("Error reading feature file");
			if(p==0)
			{
				dmin[i] = val;
				dmax[i] = val;
			}
			else
			{
				if(val<dmin[i]) dmin[i] = val;
				if(val>dmax[i]) dmax[i] = val;
			}
		}
	}
	// We reset the file to the position expected
    // rewind file
    fseek(fp, 0, SEEK_SET);
    // read in number of features
    fscanf(fp, "%d", &nFeatures);
#endif

    // load data
    for (p=0; p<nPoints; p++) {
        j=0;
        for(i=0; i<nFeatures; i++) {
            float readfloatval;
            status = fscanf(fp, "%f", &readfloatval);
            val = (scalar)readfloatval;
            if (status==EOF) Error("Error reading feature file");
#ifdef STORE_DATA_AS_INTEGER
			val = (val-dmin[i])/(dmax[i]-dmin[i]);
#endif

            if (UseFeatures[0] == 0) //when we want all the features
			{
                if(i<UseLen)
#ifdef STORE_DATA_AS_INTEGER
                    Data[p*nDims + j++] = data_int_from_scalar(val);
#else
                    Data[p*nDims + j++] = val;
#endif
            }
            else  // When we want the subset specified by the binary string UseFeatures, e.g. 111000111010101
            {
                if(i<UseLen && UseFeatures[i]=='1') 
#ifdef STORE_DATA_AS_INTEGER
                    Data[p*nDims + j++] = data_int_from_scalar(val);
#else
                    Data[p*nDims + j++] = val;
#endif
            }
        }
    }

    if(UseDistributional) //replaces if(UseFloatMasks)
    {
        // rewind file
        fseek(fpfmask, 0, SEEK_SET);

        // read in number of features
        fscanf(fpfmask, "%d", &nmaskFeatures);

        if (nFeatures != nmaskFeatures)
            Error("Error: Float Mask file and Fet file incompatible");

        // load float masks
        for (p=0; p<nPoints; p++) {
            j=0;
            for(i=0; i<nFeatures; i++)
            {
                float readfloatval;
                status = fscanf(fpfmask, "%f", &readfloatval);
                if (status==EOF) Error("Error reading fmask file");
                val = (scalar)readfloatval;

                if (UseFeatures[0] == 0)
                {
                    if(i<UseLen )
                    {
#ifdef STORE_FLOAT_MASK_AS_CHAR
						CharFloatMasks[p*nDims+j] = convert_to_char(val);
#else
                        FloatMasks[p*nDims + j] = val;
#endif
                        j++;
                    }
                }
                else  // When we want all the features
                {
                    if(i<UseLen && UseFeatures[i]=='1'  ) //To Do: implement DropLastNFeatures
                    {
#ifdef STORE_FLOAT_MASK_AS_CHAR
						CharFloatMasks[p*nDims+j] = convert_to_char(val);
#else
                        FloatMasks[p*nDims + j] = val;
#endif
                        j++;
                    }
                }
            }
        }
    }    

#ifndef COMPUTED_BINARY_MASK
    if(UseDistributional)
    {
        for(p=0; p<nPoints; p++)
            for(i=0; i<nDims; i++)
            {
#ifdef STORE_FLOAT_MASK_AS_CHAR
                if(CharFloatMasks[p*nDims+i]==(unsigned char)255) //changed so that this gives the connected component masks
#else
                if(FloatMasks[p*nDims+i]==(scalar)1) //changed so that this gives the connected component masks
#endif
                    Masks[p*nDims+i] = 1;
                else
                    Masks[p*nDims+i] = 0;
            }
    }
    else  //Case for Classical EM KlustaKwik
#else
	if(!UseDistributional)
#endif
    {
        for(p=0; p<nPoints; p++)
            for(i=0; i<nDims; i++)
                Masks[p*nDims+i] = 1;
    }

    fclose(fp);
    if(UseDistributional)
        fclose(fpfmask);

#ifndef STORE_DATA_AS_INTEGER
    // normalize data so that range is 0 to 1: This is useful in case of v. large inputs
    for(i=0; i<nDims; i++) {

        //calculate min and max
        min = HugeScore; max=-HugeScore;
        for(p=0; p<nPoints; p++) {
            val = Data[p*nDims + i];
            if (val > max) max = val;
            if (val < min) min = val;
        }

        // now normalize
        for(p=0; p<nPoints; p++) Data[p*nDims+i] = (Data[p*nDims+i] - min) / (max-min);
    }
#endif

    Output("----------------------------------------------------------\nLoaded %d data points of dimension %d.\n", (int)nPoints, (int)nDims);
    Output("MEMO: A lower score indicates a better clustering \n ");
}

// write output to .clu file - with 1 added to cluster numbers, and empties removed.
void KK::SaveOutput()
{
    integer c;
    uinteger p;
    char fname[STRLEN];
    FILE *fp;
    integer MaxClass = 0;
    vector<integer> NotEmpty(MaxPossibleClusters);
    vector<integer> NewLabel(MaxPossibleClusters);

    // find non-empty clusters
    for(c=0;c<MaxPossibleClusters;c++) NewLabel[c] = NotEmpty[c] = 0;
    for(p=0; p<BestClass.size(); p++) NotEmpty[BestClass[p]] = 1;

    // make new cluster labels so we don't have empty ones
    NewLabel[0] = 1;
    MaxClass = 1;
    for(c=1;c<MaxPossibleClusters;c++) {
        if (NotEmpty[c]) {
            MaxClass++;
            NewLabel[c] = MaxClass;
        }
    }

    // print file
    sprintf(fname, "%s.clu.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");

    fprintf(fp, "%d\n", (int)MaxClass);
    for (p=0; p<BestClass.size(); p++) fprintf(fp, "%d\n", (int)NewLabel[BestClass[p]]);

    fclose(fp);

    if(SaveCovarianceMeans)
        SaveCovMeans();
    if(SaveSorted&&UseDistributional)
        SaveSortedClu();
}

// write output to .clu file - with 1 added to cluster numbers, and empties removed.
void KK::SaveTempOutput()
{
    integer c;
    uinteger p;
    char fname[STRLEN];
    FILE *fp;
    integer MaxClass = 0;
    vector<integer> NotEmpty(MaxPossibleClusters);
   vector<integer> NewLabel(MaxPossibleClusters);
    
    // find non-empty clusters
    for(c=0;c<MaxPossibleClusters;c++) NewLabel[c] = NotEmpty[c] = 0;
   // We are merely storing the results of the current iteration, 
   //it may not be the best so far
   for(p=0; p<Class.size(); p++) NotEmpty[Class[p]] = 1;
    
    // make new cluster labels so we don't have empty ones
    NewLabel[0] = 1;
    MaxClass = 1;
    for(c=1;c<MaxPossibleClusters;c++) {
        if (NotEmpty[c]) {
            MaxClass++;
            NewLabel[c] = MaxClass;
        }
    }
    
    // print temp.clu file
   //This is the clu for the current iteration
   //This fixes the bug of having a trivial temp.clu file if there is only one iteration
    sprintf(fname, "%s.temp.clu.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");
    
    fprintf(fp, "%d\n", (int)MaxClass);
    for (p=0; p<Class.size(); p++) fprintf(fp, "%d\n", (int)NewLabel[Class[p]]);
    fclose(fp);
    
    // print besttemp.clu file
   //This is the best so far
 //   sprintf(fname, "%s.besttemp.clu.%d", FileBase, (int)ElecNo);
 //   fpb = fopen_safe(fname, "w");
    
 //   fprintf(fpb, "%d\n", (int)BestMaxClass);
 //   for (p=0; p<BestClass.size(); p++) fprintf(fpb, "%d\n", (int)BestNewLabel[BestClass[p]]);
 //   fclose(fpb);
   
   
    if(SaveCovarianceMeans)
        SaveCovMeans();
    if(SaveSorted&&UseDistributional)
        SaveSortedClu();
}

void KK::SaveCovMeans()
{
    char fname[STRLEN];
    FILE *fp;
    //// print covariance to file
    //sprintf(fname, "%s.cov.%d", FileBase, (int)ElecNo);
    //fp = fopen_safe(fname, "w");
    //for (integer cc=0; cc<nClustersAlive; cc++)
    //{
    //    integer c = AliveIndex[cc];
    //    for(integer i=0; i<nDims; i++)
    //    {
    //        for(integer j=0; j<nDims; j++)
    //        {
				//// TODO: update Cov output for distributional
				//if (!UseDistributional)
				//	fprintf(fp, SCALARFMT " ", Cov[c*nDims2+i*nDims+j]);
    //        }
    //        fprintf(fp, "\n");
    //    }
    //    fprintf(fp, "\n");
    //}
    //fclose(fp);
    // print mean to file
    sprintf(fname, "%s.mean.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");
    for (integer cc=0; cc<nClustersAlive; cc++)
    {
        integer c = AliveIndex[cc];
        for(integer i=0; i<nDims; i++)
        {
            fprintf(fp, SCALARFMT " ", Mean[c*nDims+i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// Saves sorted.fet and sorted.mask file
void KK::SaveSortedData()
{
    char fname[STRLEN];
    FILE *fp;
    // sorted.fet file
    sprintf(fname, "%s.sorted.fet.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");
    fprintf(fp, "%d\n", (int)nDims);
    for(integer q=0; q<nPoints; q++)
    {
        integer p = SortedIndices[q];
        for(integer i=0; i<nDims; i++)
            fprintf(fp, SCALARFMT " ", GetData(p, i));
        fprintf(fp, "\n");
    }
    fclose(fp);
    // sorted.mask file
    sprintf(fname, "%s.sorted.mask.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");
    fprintf(fp, "%d\n", (int)nDims);
    for(integer q=0; q<nPoints; q++)
    {
        integer p = SortedIndices[q];
        for(integer i=0; i<nDims; i++)
            fprintf(fp, "%d ", (int)GetMasks(p*nDims+i));
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// Save sorted.clu file (see SaveOutput for explanation)
void KK::SaveSortedClu()
{
    char fname[STRLEN];
    FILE *fp;
    vector<integer> NotEmpty(MaxPossibleClusters);
    vector<integer> NewLabel(MaxPossibleClusters);
    for(integer c=0; c<MaxPossibleClusters; c++)
        NewLabel[c] = NotEmpty[c] = 0;
    for(integer q=0; q<nPoints; q++)
        NotEmpty[Class[SortedIndices[q]]] = 1;
    NewLabel[0] = 1;
    integer MaxClass = 1;
    for(integer c=1; c<MaxPossibleClusters; c++)
        if(NotEmpty[c])
            NewLabel[c] = ++MaxClass;
    sprintf(fname, "%s.sorted.clu.%d", FileBase, (int)ElecNo);
    fp = fopen_safe(fname, "w");
    fprintf(fp, "%d\n", (int)MaxClass);
    for(integer q=0; q<nPoints; q++)
        fprintf(fp, "%d\n", (int)NewLabel[Class[SortedIndices[q]]]);
    fclose(fp);
}
