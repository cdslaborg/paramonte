!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This is main entry to the tests of the ParaMonte kernel library.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program main

use pm_test, only: setSummary

call random_seed()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

block; use test_pm_arrayCenter; call setTest(); end block
block; use test_pm_arrayChange; call setTest(); end block
block; use test_pm_arrayChoice; call setTest(); end block
block; use test_pm_arrayCompact; call setTest(); end block
block; use test_pm_arrayComplement; call setTest(); end block
block; use test_pm_arrayCompareLex; call setTest(); end block
block; use test_pm_arrayCopy; call setTest(); end block
block; use test_pm_arrayCopy; call setTest(); end block
block; use test_pm_arrayFill; call setTest(); end block
block; use test_pm_arraySearch; call setTest(); end block
block; use test_pm_arrayFind; call setTest(); end block
block; use test_pm_arrayInit; call setTest(); end block
block; use test_pm_arrayInsert; call setTest(); end block
block; use test_pm_arrayMerge; call setTest(); end block
block; use test_pm_arrayPad; call setTest(); end block
block; use test_pm_arrayRange; call setTest(); end block
block; use test_pm_arrayRank; call setTest(); end block
block; use test_pm_arrayRefill; call setTest(); end block
block; use test_pm_arrayRebill; call setTest(); end block
block; use test_pm_arrayRebind; call setTest(); end block
block; use test_pm_arrayRemap; call setTest(); end block
block; use test_pm_arrayRemove; call setTest(); end block
block; use test_pm_arrayReplace; call setTest(); end block
block; use test_pm_arrayResize; call setTest(); end block
block; use test_pm_arrayReverse; call setTest(); end block
block; use test_pm_arraySpace; call setTest(); end block
block; use test_pm_arraySelect; call setTest(); end block
block; use test_pm_arrayShuffle; call setTest(); end block
block; use test_pm_arraySort; call setTest(); end block
block; use test_pm_arraySplit; call setTest(); end block
block; use test_pm_arrayStrip; call setTest(); end block
block; use test_pm_arrayUnique; call setTest(); end block
block; use test_pm_arrayVerbose; call setTest(); end block
block; use test_pm_bench; call setTest(); end block
block; use test_pm_complexAbs; call setTest(); end block
block; use test_pm_complexDiv; call setTest(); end block
block; use test_pm_complexCompareAll; call setTest(); end block
block; use test_pm_complexCompareAny; call setTest(); end block
block; use test_pm_complexCompareLex; call setTest(); end block
block; use test_pm_cosmicRate; call setTest(); end block
block; use test_pm_cosmology; call setTest(); end block
block; use test_pm_dateTime; call setTest(); end block
block; use test_pm_distanceEuclid; call setTest(); end block
block; use test_pm_distBern; call setTest(); end block
block; use test_pm_distExp; call setTest(); end block
block; use test_pm_distGamma; call setTest(); end block
block; use test_pm_distGenExpGamma; call setTest(); end block
block; use test_pm_distPareto; call setTest(); end block
block; use test_pm_distPower; call setTest(); end block
block; use test_pm_distPiwiPoweto; call setTest(); end block
block; use test_pm_except; call setTest(); end block
block; use test_pm_mathCompare; call setTest(); end block
block; use test_pm_mathCumPropExp; call setTest(); end block
block; use test_pm_mathCumSum; call setTest(); end block
block; use test_pm_mathFactorial; call setTest(); end block
block; use test_pm_mathFactoring; call setTest(); end block
block; use test_pm_mathExp; call setTest(); end block
block; use test_pm_mathRoot; call setTest(); end block
block; use test_pm_matrixChol; call setTest(); end block
block; use test_pm_matrixMulAdd; call setTest(); end block
block; use test_pm_matrixMulTri; call setTest(); end block
block; use test_pm_sampleCCF; call setTest(); end block
block; use test_pm_sampleCor; call setTest(); end block
block; use test_pm_sampleCov; call setTest(); end block
block; use test_pm_sampleMean; call setTest(); end block
block; use test_pm_sampleShift; call setTest(); end block
block; use test_pm_sampleVar; call setTest(); end block
block; use test_pm_timer; call setTest(); end block

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!block; use test_pm_matrixInit; call setTest(); end block

!!block; use test_pm_distBand; call setTest(); end block
!!block; use test_pm_batse; call setTest(); end block
!!block; use test_pm_except; call setTest(); end block
!!block; use test_pm_val2complex; call setTest(); end block
!!block; use test_pm_sampleCrossCorr; call setTest(); end block
!!block; use test_pm_distanceMahal; call setTest(); end block
!!block; use test_pm_distMultiNorm; call setTest(); end block
!!block; use test_pm_distMultiNorm; call setTest(); end block
!!block; use test_pm_distUnifEll; call setTest(); end block
!!block; use test_pm_distMultiSkewNorm; call setTest(); end block
!!block; use test_pm_distNorm; call setTest(); end block

!!block; use test_pm_distUnif; call setTest(); end block
!!block; use test_pm_distUnif; call setTest(); end block
!!block; use test_pm_domainCube; call setTest(); end block
!!block; use test_pm_domainBall; call setTest(); end block
!!block; use test_pm_err; call setTest(); end block
!!block; use test_pm_io; call setTest(); end block
!!block; use Test_FileContents_pmod; call setTest(); end block
!!block; use test_pm_io; call setTest(); end block
!!block; use test_pm_distGeomCyclic; call setTest(); end block
!!block; use test_pm_knn; call setTest(); end block
!!block; use test_pm_ellipsoid; call setTest(); end block
!!block; use test_pm_histogram; call setTest(); end block
!!block; use test_pm_statest; call setTest(); end block
!!block; use test_pm_val2int; call setTest(); end block
!!!block; use test_pm_polation; call setTest(); end block
!!block; use test_pm_quadRomb; call setTest(); end block
!!block; use test_pm_polation; call setTest(); end block
!!block; use test_pm_clusKmeans; call setTest(); end block
!!block; use test_pm_logicalCompare; call setTest(); end block
!!block; use test_pm_val2logical; call setTest(); end block
!
!!block; use test_pm_math; call setTest(); end block
!!block; use test_pm_mathGamma; call setTest(); end block
!!block; use test_pm_mathLogAddExp; call setTest(); end block
!!block; use test_pm_mathLogSumExp; call setTest(); end block
!!block; use test_pm_mathMinMax; call setTest(); end block
!
!block; use test_pm_matrix; call setTest(); end block
!block; use test_pm_matrixDet; call setTest(); end block
!block; use test_pm_matrixInit; call setTest(); end block
!block; use test_pm_matrixCopy; call setTest(); end block
!block; use test_pm_matrixTrans; call setTest(); end block
!
!block; use test_pm_optimization; call setTest(); end block
!block; use test_pm_option; call setTest(); end block
!block; use test_pm_partition; call setTest(); end block

!!!block; use test_pm_partitionRecursiveMaxDen; call setTest(); end block
!!!block; use test_pm_partitionRecursiveMinVol; call setTest(); end block
!!!block; use Test_PartitionBenchm_pmod; call setTest(); end block
!!!block; use test_pm_partitionMaxDen; call setTest(); end block

!block; use test_pm_parallelism; call setTest(); end block
!block; use test_pm_sysPath; call setTest(); end block
!!block; use test_pm_processPoisson; call setTest(); end block
!block; use test_pm_sampleQuan; call setTest(); end block
!block; use test_pm_randomSeed; call setTest(); end block
!block; use test_pm_val2real; call setTest(); end block
!block; use Test_SampleCovMat_pmod; call setTest(); end block
!block; use test_pm_sampleECDF; call setTest(); end block
!block; use test_pm_sampleShift; call setTest(); end block
!block; use Test_SampleVariance_pmod; call setTest(); end block
!block; use Test_Set_pmod; call setTest(); end block
!block; use test_pm_statistics; call setTest(); end block
!!block; use test_pm_mathGamma; call setTest(); end block
!block; use test_pm_strASCII; call setTest(); end block
!block; use test_pm_container; call setTest(); end block
!block; use test_pm_val2str; call setTest(); end block
!block; use test_pm_str; call setTest(); end block
!block; use test_pm_sysShell; call setTest(); end block
!block; use test_pm_timer; call setTest(); end block
!block; use test_pm_timer; call setTest(); end block
!block; use test_pm_tranGaus; call setTest(); end block

#if SAMPLER_TEST_ENABLED && 1
!block; use test_pm_paraDRAM; call setTest(); end block
!block; use test_pm_paraDISE; call setTest(); end block
!block; use test_pm_paraNest; call setTest(); end block
!block; use test_pm_paraDRAM_RefinedChain; call setTest(); end block
!block; use test_pm_paraDISE_RefinedChain; call setTest(); end block
!block; use test_pm_paraDRAM_ChainFileContents; call setTest(); end block
!block; use test_pm_paraDISE_ChainFileContents; call setTest(); end block
!!block; use test_pm_paraNest_ChainFileContents; call setTest(); end block
#endif

call setSummary()

!block
!    use pm_except, only: BELL
!    write(*,*) BELL
!end block

end program main