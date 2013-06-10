KlustaKwik
==========

1) Introduction
------------------------
------------------------

KlustaKwik is a program for cluster analysis by fitting a mixture of Gaussians. It was designed for the specific problem of spike sorting of multi-electrode arrays, but can be used for any application. Technically, KlustaKwik works by implementing a hard EM algorithm with unconstrained covariance matrices. It also uses a number of tricks to greatly speed up execution. It has recently been updated to deal with high-dimensional data as required for spike sorting from large, dense electrode arrays. An older version of the code, which cannot deal with high-density arrays, can be found [here](http://sourceforge.net/projects/klustakwik/). But don't use that version, use this one: even for plain tetrode data, the new version is fully backward compatible, and is also much faster.

For spike sorting, KlustaKwik is designed to be used in conjunction with [SpikeDetekt](http://klusta-team.github.io/spikedetekt) for spike detection, and [KlustaViewa](http://klusta-team.github.io/klustaviewa) for manual verification and adjustment of clustering results.

2) Installation
---------------------
---------------------

KlustaKwik is written in C++, and provided as source code. A makefile is provided. It is a plain-text program that runs from the command line, so compilation should be straightforward. 

On Linux and Mac OS X simply unzip the source folder, open a command line terminal in the unzipped folder and type:

	make
	
This will create the executable: KlustaKwik.

3) Usage
---------------------
---------------------

KlustaKwik 3.0 is backward compatible with previous versions. To use it in "classic" mode, just run the same command you would have for version 2.x (as documented on the sourceforge page linked above). It will produce the same results, but should run about 10 times faster. In this mode, KlustaKwik takes a single input file (mydata.fet.n) containing feature vectors, and produces an output file (mydata.clu.n) containing cluster numbers. It also produces a log file (mydata.klg.n).

To deal with high channel count arrays, KlustaKwik should be run in "masked" mode. In this mode, the algorithm takes advantage of the fact that any spike occurs on only a subset of channels, with the remainder containing only multi-unit noise. The information of which channels are relevant for each spike is encoded in an additional .fmask file which is output from SpikeDetekt, along with the usual .fet features file.  *Unmasked* channels are channels on which spiking activity has been found to occur by the program SpikeDetekt, whereas *masked* channels contain only noise. The .fmask file is a text file, every line of which is a vector of length the number of features, in which **1** denotes **unmasked** and **0** denotes **masked**, and values between 0 and 1 indicate partial masking.

The suffix .n allows you to keep track of data recorded on a probe with mutliple shanks. Thus, for shank n, the relevant files would be the following output from SpikeDetekt:

    mydata.fet.n
    mydata.fmask.n


4) Command line input
----------------------
----------------------

KlustaKwik runs from the command line, and takes a large number of options. We plan to rationalize these, but in the meantime you have to run fairly long command strings to run it in masked mode.

The first two arguments are the filebase and shank number, after which come optional parameters. For example If you wanted to cluster the 4th shank from a file called "recording" in masked mode, you would run something like this:

    [yourterminal]$./KlustaKwik recording 4 -UseDistributional 1 -UseMaskedInitialConditions 1 -AssignToFirstClosestMask 1 -MaxPossibleClusters 500 -MinClusters 130 -MaxClusters 130 -PenaltyK 1 -PenaltyKLogN 0 -UseFeatures 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110

You will probably want to write a script to generate such these commands. An example Python script could read as follows:

    import os, sys

    filebase = 'recording'
    shank_num = '4'
    num_features = 97
    #Number of features (including time)
    # Replace './KlustaKwik' with the path on your system pointing to the executible KlustaKwik
    
    os.system(
        	'./KlustaKwik'
		    ' '+filebase+' '+shank_num+' -UseFeatures '+'1'*(num_features-1)+'0'+' '
    		'-MinClusters 200 '
            '-MaxClusters 200 '
    		'-MaxPossibleClusters 500 '
    		'-PenaltyK 1.0 '
            '-PenaltyKLogN 0.0 '
    		'-PriorPoint 1 '
    		'-Debug 0 '
    		'-RandomSeed 654 '
    		'-SaveSorted 0 '
    		'-UseMaskedInitialConditions 1 '
            '-AssignToFirstClosestMask	1'
    		'-UseDistributional 1 '
    		)

You would then edit the parameters filebase, shank_num, and num_features as required. 

5) Input Files
--------------
--------------

The first line of the .fet file is a single number specifying the number of features. After this, each line contains a vector giving the feature values of one spike. For example, for 32 channels, 3 principal components per channel, and a final feature representing time, the feature file would start off like:

    97
    -43206 -6407 652 5181 3675 11860 -30066 -145 -3720 -33491 -15892 5840 -57344 11025 -16524 -8815 -4046 -157 -30493 -7224 507 -676 9700 -6809 -39397 -211 -5468 -22511 -352 -186 -1127 -10679 9435 -25851 -22021 4238 -26672 -17418 11085 10035 -2864 9975 -20875 -29624 13367 -14849 -8633 4408 -19751 -24882 8184 -28165 4853 2136 -21114 -10823 3004 -8716 -4686 4793 11529 -9505 8930 13056 5133 3442 -1085 -4045 4965 -63934 -3827 -6996 -46672 4384 -7775 -29050 -549 -9260 -24985 20865 -14532 -17785 -6168 862 -22615 1593 -5319 -972 -1660 10063 -106 1600 5378 1963 -4880 7357 2
    -43496 -2123 1697 7471 2085 15302 -30109 2679 -2230 -34578 -13000 6082 -56265 17120 -14871 -9465 -2789 -1098 -31967 -3309 -255 -74 10289 -6540 -39518 3964 -4349 -22094 1344 1581 -1050 -11070 9237 -28378 -19400 2499 -27143 -15752 11980 10508 -3931 10804 -22956 -28054 11357 -15002 -8109 5675 -22333 -22982 6313 -26356 7494 4534 -21650 -8983 2692 -8909 -4412 5892 10953 -10840 8443 14319 3268 5287 -236 -4808 6393 -64397 2215 -5170 -45979 8712 -5203 -29763 2625 -9344 -22857 23174 -11572 -17982 -5042 1361 -22808 4433 -5423 1133 -3507 14082 1741 292 7829 2157 -5736 8542 2
    -33069 18324 -563 -108726 -20512 -11860 -30497 10828 20051 -31983 -2309 24404 -11377 22760 14372 -27932 -2754 7252 -9373 43163 -5022 -48048 22884 -11330 -15395 36739 28855 -14419 13872 -976 -45777 305 5922 -4145 9977 5397 -9010 -3391 -13322 7835 278 30982 -47770 -15958 -364 -22770 1765 3022 -1997 8759 13377 -42757 -19253 27619 6931 -4810 1240 -14390 -39503 -9258 2612 -25492 -7973 -16064 -13039 -3647 -17640 -11906 11890 2259 -8161 -12922 -23660 -28462 7323 -18757 -25355 -6518 -15544 -26716 12220 3591 -9433 -32399 -5565 -36616 -7977 18442 -10641 -28989 -14625 -404 -1207 6840 -17422 -3009 41
    2903 -6480 -8171 -6906 -14576 10433 29882 -4122 -7722 6632 29771 9732 -5294 -3917 50886 -7037 2175 4469 11817 -4083 -7448 -7708 3939 3240 -14890 2721 4889 -11259 20546 11612 -36244 -5598 27674 9815 -4252 -16080 -13403 6912 11855 8481 11759 -24013 -41589 -36261 11447 -32228 2745 -16386 -9704 4388 -20931 3564 32431 -9779 -5692 26247 -11884 -1589 37336 12236 3088 32036 -18553 14636 18875 1886 22541 10968 -14971 8588 23946 -23505 -26785 70034 -29746 -28573 50845 -24927 -166980 114185 -57570 -83298 67056 -39670 -103696 90495 -36092 -52810 31849 -5378 -6437 3499 -34861 -11266 -8738 24618 84

The .fmask file has a similar format, but with each line giving a vector of floats between 0 and 1. The .fmask file also has the first line specifying the number of features. 

6) Parameters
-------------------
-------------------

There are a large number of parameters that control how KlustaKwik works. The defaults, which cause it to run in classic mode, are:

    FileBase    electrode
    ElecNo	1
    MinClusters	20
    MaxClusters	30
    MaxPossibleClusters	100
    nStarts	1
    RandomSeed	1
    Debug	0
    Verbose	1
    UseFeatures	11111111111100001
    DistDump	0
    DistThresh	6.907755
    FullStepEvery	20
    ChangedThresh	0.050000
    Log	1
    Screen	1
    MaxIter	500
    StartCluFile	
    SplitEvery	40
    PenaltyK	0.000000
    PenaltyKLogN	1.000000
    Subset	1
    PriorPoint	1
    SaveSorted	0
    SaveCovarianceMeans	0
    UseMaskedInitialConditions	0
    AssignToFirstClosestMask	0
    UseDistributional	0
    help	0


The most important parameters are: 

+ **UseFeatures** 

After this comes a string with 1's for features you want to use, and 0's for features you don't want to use. In classic mode, you use this option to take out bad channels. In masked mode, you should instead take bad channels out from the .probe file. The default is a historical string for tetrodes - you should always override it.

+ **Penalties**

KlustaKwik uses penalties to reduce the number of clusters fit. The parameters PenaltyK and PenaltyKLogN are given positive values. The higher the values, the fewer clusters you obtain. Higher penalties discourage cluster splitting. PenaltyKLogN also increases penalty when there are more points.

-PenaltyK 0 -PenaltyKLogN 1 is the default, corresponding to the "Bayesian Information Criterion". This is recommended for classic mode.

-PenaltyK 1 -PenaltyKLogN  corresponds to "Akaike's Information Criterion". This is recommended for large probes and masked mode.

+ **UseDistributional** (*default* 0)

To use KlustaKwik in "masked" mode, set this to 1.
This enables the use of the new `masked Expectation-Maximization' algorithm. To use the new algorithm, it is recommended that you change most of the defaults.

+ **MinClusters**, **MaxClusters**

In classic mode, KlustaKwik starts from random cluster assignments, running a new random start for every integer between MinClusters and MaxClusters. 

+ **UseMaskedInitialConditions**, **AssignToFirstClosestMask** 

In masked mode, a better approach is to start with clusters that are derived from the mask vectors. Set both of these options to 1 to enable masked initializations. This gives better performance and runs faster. Also set **MinClusters**= **MaxClusters** to the number of distinct clusters it will start from, and it will assign each spike to the nearest unique mask according to Hamming distance.

**SplitEvery** (*default* 40)
is an integer which is the number of iterations after which KlustaKwik attempts to split existing clusters.
When using masked initializations, to save time due to excessive splitting, set **SplitEvery** to a large number, close
to the number of distinct masks or the number of chosen starting masks. 


7) Full Glossary of Parameters
-------------------------
-------------------------
**UseDistributional (default 0)** - Set this to 1 to use in "masked" mode.

**Filebase** - Name of your .fet and .mask file, e.g. if your feature file is called mydata.fet.1, then Filebase is *mydata*.

**ElecNo** - Shank number of your probe.

**MinClusters n** (default 20)  - The minimum number of starting clusters. It will give a random initial assignment when the **-UseMaskedInitialConditions** is not being used; otherwise it will assign according to mask. The intial assignment will have no less than n clusters.  The final number may be different, since clusters can be split or deleted during the course of the algorithm

**MaxClusters n** (default 30) - The maximum number of starting clusters. It will give a random initial assignment when the **-UseMaskedInitialConditions** is not being used; otherwise it assignment will be determined by masks. Note: It ought to be set lower than **MaxPossibleClusters**.
 
**nStarts n** (default 1) -  The algorithm will be started n times for each inital cluster count between MinClusters and MaxClusters.

**MaxPossibleClusters n** (default 100) - The largest permitted number of clusters, so cluster splitting can produce no more than n clusters. Note: It ought to be set higher than **MaxClusters**.

**UseMaskedInitialConditions	(default 0)** - Initialises using distinct derived binary masks. Use together with **AssignToFirstClosestMask** below. See previous section for explanation.

**AssignToFirstClosestMask	(default 0)** - If starting with a number of clusters fewer than the number of distinct derived binary masks, it will assign the rest of the points to the cluster with the nearest mask.

**help** - Prints a short message and then the default parameter values.

**RandomSeed n**    (default 1) Specifies a seed for the random number generator.

**StartCluFile** STRING   (default " ") Initializes according to the specified cluster file.  If it can't find a better cluster assignment, it will output this.

**FullStepEvery** n (default 10) All log-likelihoods are recalculated every n steps (see DistThresh).

**SplitEvery** n    (default 40) Test to see if any clusters should be split every n steps. 0 completely suppresses splitting altogether.

**MaxIter** n       (default 500) Maximum number of iterations. ie. it won't try more than n iterations from any starting point.


Various debugging options for developers:
-----------------------------------------

**Log**             (default 1) Produces .klg log file (default is yes, to switch off do -Log 0).

**Screen**          (default 1) Produces parameters and progress information on the console. Set to 0 to suppress output in batches.

**SaveSorted** (default 0)	Saves a .clu file with masks sorted lexicographically.
    
**SaveCovarianceMeans**	(default 0) Saves means and covariance matrices. Stops computation at each iteration. Manual input required for continuation.

**Debug**           (default 0) Miscellaneous debugging information (not recommended).

**DistDump**        (default 0) Outputs a ridiculous amount of debugging information (definately not recommended).

**PriorPoint** (default 1) Helps normalize covaraince matrices.

**DistThresh d**     (default 6.907755) Time-saving parameter. If a point has log likelihood more than d worse for a given class than for the best class, the log likelihood for that class is not recalculated.  This saves an awful lot of time.

**ChangedThresh f**  (default 0.05) All log-likelihoods are recalculated if the fraction of instances changing class exceeds f (see DistThresh).
 
8) Output Files
---------------
---------------

Regardless of the parameters chosen to run KlustaKwik, the two following plain text files will always be produced.

	mydata.clu.n (main output file)
	mydata.klg.n (a log file)

The .clu file will be a file with with S+1 lines, where S is the number of spikes detected. The first line will give the total number of clusters found. Each subsequent line will give the cluster label corresponding to each spike.e.g. if the start of your .clu file looks like this:

	101
	15
	3
	4
	79
	51
	29
	4
        .
        .

This will mean there are 101 clusters in total. That the 1st is in cluster 15, the 2nd in cluster 3 and the 3rd and 7th spike both in cluster 4, etc. until the end of the file.

If you are uncertain whether or not you have run KlustaKwik having inputted the correct parameters simply open the .klg file (it is produced as soon as you run the program) and all parameter values will be listed. 


	
	




