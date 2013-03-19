KlustaKwik
==========

Documentation for Masked KlustaKwik
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

0) Introduction
------------------------
------------------------

The new masked version of KlustaKwik, currently KlustaKwik3.0.2, differs from previous versions in its use of a maskfile (either binary in the form of a .mask file  or floats between 0 and 1, .fmask file) in addition to the usual features file (.fet file). 

It is designed to be used in conjunction with SpikeDetekt for clustering spike waveforms recorded on large dense probes. 

*Unmasked* channels are channels on which spiking activity has been found to occur by the program SpikeDetekt, whereas *masked* channels contain only noise.The .mask file or the .fmask file is a text file, every line of which is a vector giving the positions of the unmasked channels. In the .mask file, **1** denotes *unmasked* and **0** denotes *masked*. The .fmask file in contrast permits values between 0 and 1.

1) Parameters
-------------------
-------------------

The current release of masked KlustaKwik has an enormous range of parameters which can be adjusted according to the user's needs. Many unnecessary parameters will soon be phased out, but in the interim to help beta-testers, we will describe a few of the relevant ones.

    Usage: MaskedKlustaKwik FileBase ElecNo [Arguments]

    Arguments (with default values): 

    FileBase  electrode
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
    UseMaskedEStep	0
    UseMaskedMStep	0
    UseClusterPenalty	0
    SaveCovarianceMeans	0
    UseMaskedInitialConditions	0
    AssignToFirstClosestMask	0
    UseDistributional	0
    UseDistributionalEStep	1
    UseFloatMasks	0
    help	0

The above defaults cause KlustaKwik to run exactly as previous versions, but 10 times faster. The only difference is that the parameter PenaltyMix has been replaced with two parameters, PenaltyK and PenaltyKLogN.

+ **Penalties**

 

PenaltyMix 1 corresponds to (PenaltyK 0, PenaltyKLogN 1)(Bayesian Information Criterion)

PenaltyMix 0 corresponds to (PenaltyK 1, PenaltyKLogN 0) (Akaike Information Criterion)

The parameters PenaltyK and PenaltyKLogN can be given positive values. The higher the values, the fewer clusters you obtain. Higher penalties discourage splitting.

AIC is recommended for larger probes.

+ **UseDistributional, UseDistributionalEStep and UseClusterPenalty**

These enables the use of the new `distributional Expectation-Maximization' algorithm. Set these all to 1. 

+ **UseFloatMasks**

It has been observed that using soft masks of the form. e.g.:

0 0 0 0.3 0.3 0.3 1 1 1 1 1 1 0.4 0.4 0.4 0 0 0 0 0 0 

leads to improved clusterings than using binary masks:

0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 

Set this to 1 to enable masks to take floating point values. 

+ **UseMaskedInitialConditions**, **AssignToFirstClosestMask** and **SplitEvery**

The quality of clustering obtained from EM-type algorithms is sensitive to the initialization. A good way to initialize that works better than random initializations is to start with clusters which depend on unique masks or clusters containing similar masks. Set both of these options to 1 to enable masked initializations. At the same time set **MinClusters**= **MaxClusters** to the number of distinct masks (also set MaxPossibleClusters to a value greater than this). 

(If you do not know the number of distinct masks, run KlustaKwik, with arbitrary values of MinClusters and MaxClusters, read the first line giving the number of distinct masks, and then prematurely terminate the program).

**SplitEvery** is an integer which is the number of iterations after which KlustaKwik attempts to split existing clusters. When using masked initializations, to save time due to excessive splitting, set **SplitEvery** to a large number, close to the number of distinct masks. 

+ **UseFeatures**

Make sure to include a **1** for every feature you would like to include and a **0** for every feature you want to leave out (i.e. features corresponding to bad channels that you don't want).

+ ***Warning: Ignore UseMaskedEStep and UseMaskedMStep***

They should be set to 0 (as is the default). They are obsolete.

2) Command line input
----------------------
----------------------

A typical command to run the masked version of KlustaKwik therefore looks as follows in a linux terminal:

    [yourterminal]$./KlustaKwik yourfetfilename shanknumber -UseDistributional 1 -UseDistributionalEStep 1 -UseMaskedMStep 0 -UseMaskedEStep 0 -FullStepEvery 1 -SplitEvery 40 -UseMaskedInitialConditions 1 -AssignToFirstClosestMask 1 -UseClusterPenalty 1 -UseFloatMasks 0 -RandomSeed 123 -MaxIter 10000 -MaxPossibleClusters 500 -MinClusters 130 -MaxClusters 130 -PenaltyK 1 -PenaltyKLogN 0 -UseFeatures 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110
    
e.g. if your .fet file is called **recording.fet.4** (and your other files are **recording.mask.4**, **recording.fmask.4**)  for the fourth shank, then the command looks like:

    [yourterminal]$./KlustaKwik recording 4 -UseDistributional 1 -UseDistributionalEStep 1 -UseMaskedMStep 0 -UseMaskedEStep 0 -FullStepEvery 1 -SplitEvery 40 -UseMaskedInitialConditions 1 -AssignToFirstClosestMask 1 -UseClusterPenalty 1 -UseFloatMasks 0 -RandomSeed 123 -MaxIter 10000 -MaxPossibleClusters 500 -MinClusters 130 -MaxClusters 130 -PenaltyK 1 -PenaltyKLogN 0 -UseFeatures 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110

You may consider writing a script to generate such a complicated command.

**We apologize for the current somewhat complicated set-up. Everything will be simplified once beta testing has been completed.**
