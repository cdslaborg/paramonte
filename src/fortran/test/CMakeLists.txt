####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################

# test source files

set(ParaMonteTest_SOURCE_FILES
    ${paramonte_src_fortran_test_dir}/main.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCenter.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayChange.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayChoice.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCompact.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCompareLex.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayComplement.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCopy.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayFill.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySearch.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayFind.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayInit.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayInsert.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayMerge.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayPad.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRange.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRank.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRebill.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRebind.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRefill.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRemap.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRemove.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayReplace.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayResize.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayReverse.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySelect.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayShuffle.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySort.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySpace.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySplit.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayStrip.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayUnique.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayVerbose.F90

    ${paramonte_src_fortran_test_dir}/test_pm_bench.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexAbs.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexDiv.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAll.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAny.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareLex.F90
    ${paramonte_src_fortran_test_dir}/test_pm_cosmicRate.F90
    ${paramonte_src_fortran_test_dir}/test_pm_cosmology.F90

    ${paramonte_src_fortran_test_dir}/test_pm_dateTime.F90
   #${paramonte_src_fortran_test_dir}/test_pm_io.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distBern.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distExp.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distGamma.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distGenExpGamma.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPiwiPoweto.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPower.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPareto.F90
    ${paramonte_src_fortran_test_dir}/test_pm_except.F90

    ${paramonte_src_fortran_test_dir}/test_pm_mathCompare.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathCumSum.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathCumPropExp.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathExp.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathFactorial.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathFactoring.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathRoot.F90

    ${paramonte_src_fortran_test_dir}/test_pm_matrixMulAdd.F90
    ${paramonte_src_fortran_test_dir}/test_pm_matrixMulTri.F90

   #${paramonte_src_fortran_test_dir}/test_pm_matrixInit.F90
   #${paramonte_src_fortran_test_dir}/test_pm_quadRomb.F90
#   ${paramonte_src_fortran_test_dir}/test_pm_val2complex.F90

   #${paramonte_src_fortran_test_dir}/test_pm_matrixInit@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCenter@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayChange@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayChoice@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCompact@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCompareLex@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayComplement@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayCopy@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayFill@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySearch@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayFind@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayInit@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayInsert@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayMerge@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayPad@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRange@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRank@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRebill@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRebind@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRefill@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRemap@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayRemove@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayReplace@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayResize@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayReverse@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySelect@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayShuffle@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySort@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySpace@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arraySplit@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayStrip@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayUnique@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_arrayVerbose@routines.F90

    ${paramonte_src_fortran_test_dir}/test_pm_complexAbs@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexDiv@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAll@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAny@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_complexCompareLex@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_cosmicRate@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_cosmology@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_dateTime@routines.F90

   #${paramonte_src_fortran_test_dir}/test_pm_io@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distanceEuclid.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distanceEuclid@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distBern@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distExp@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distGamma@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distGenExpGamma@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPiwiPoweto@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPareto@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_distPower@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_except@routines.F90

    ${paramonte_src_fortran_test_dir}/test_pm_mathCompare@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathCumSum@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathCumPropExp@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathExp@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathFactorial@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathFactoring@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_mathRoot@routines.F90

    ${paramonte_src_fortran_test_dir}/test_pm_matrixMulAdd@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_matrixMulTri@routines.F90

#    ${paramonte_src_fortran_test_dir}/test_pm_val2complex@routines.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_distanceMahal.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distanceMahal@routines.F90
#
#

#    ${paramonte_src_fortran_test_dir}/test_pm_distMultiNorm.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distMultiNorm.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distMultiNorm@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distUnifEll.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distMultiSkewNorm.F90
#   #${paramonte_src_fortran_test_dir}/Test_DistMultiUnif_pmod.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_distNorm.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_distUnif.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_domainCube.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_domainBall.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_domainBall@routines.F90
#
#
#    ${paramonte_src_fortran_test_dir}/test_pm_batse.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distBand.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_ellipsoid.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_err.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_io.F90
#    ${paramonte_src_fortran_test_dir}/Test_FileContents_pmod.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_distGeomCyclic.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_knn.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_histogram.F90
#    ${paramonte_src_fortran_test_dir}/pm_statest.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2int.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2int@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_polation.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_polation.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_clusKmeans.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_logicalCompare.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_logicalCompare@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2logical.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2logical@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_processPoisson.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_math.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_math@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathGamma.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathGamma@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathLogAddExp.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathLogAddExp@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathLogSumExp.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathLogSumExp@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathMinMax.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_mathMinMax@routines.F90
#
    ${paramonte_src_fortran_test_dir}/test_pm_matrixChol.F90
    ${paramonte_src_fortran_test_dir}/test_pm_matrixChol@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrixDet.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrixDet@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrixTrans.F90
#   #${paramonte_src_fortran_test_dir}/test_pm_matrixTrans@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrix.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_matrix@routines.F90
#
#    ${paramonte_src_fortran_test_dir}/test_pm_partition.F90
#    #${paramonte_src_fortran_test_dir}/Test_PartitionKmeans_pmod.F90
#    #${paramonte_src_fortran_test_dir}/Test_PartitionVarDen_pmod.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_partitionMaxDen.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_partitionRecursiveMaxDen.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_partitionRecursiveMinVol.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_misc.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_optimization.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_option.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_option@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_parallelism.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_paraDRAM_ChainFileContents.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_paraDISE_ChainFileContents.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_paraNest_ChainFileContents.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_paraDRAM_RefinedChain.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_paraDISE_RefinedChain.F90
#    #${paramonte_src_fortran_test_dir}/test_pm_paraNest_RefinedChain.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_sysPath.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_sampleQuan.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_randomSeed.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2real.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2real@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_sampleECDF.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_sampleECDF@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleMean.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleMean@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleShift.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleShift@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCCF.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCCF@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCor.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCor@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCov.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleCov@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleVar.F90
    ${paramonte_src_fortran_test_dir}/test_pm_sampleVar@routines.F90
    ${paramonte_src_fortran_test_dir}/test_pm_timer.F90
    ${paramonte_src_fortran_test_dir}/test_pm_timer@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_statistics.F90
#    ${paramonte_src_fortran_test_dir}/Test_Set_pmod.F90
#    ${paramonte_src_fortran_test_dir}/Test_Set_pmod@Subtract_smod.F90
#   #${paramonte_src_fortran_test_dir}/test_pm_mathGamma.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_str.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_str@routines.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_strASCII.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_container.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2str.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_val2str@routines.F90
#    ${paramonte_src_fortran_test_dir}/Test_pm_sysShell.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_timer.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_timer.F90
#    ${paramonte_src_fortran_test_dir}/test_pm_tranGaus.F90
#    #${paramonte_src_fortran_test_dir}/Test_EconomicsToolbox_pmod.F90
)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayCenter@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayCenter@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayChange@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayChange@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayChoice@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayChoice@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayCompact@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayCompact@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayCompareLex@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayCompareLex@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayComplement@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayComplement@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayCopy@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayCopy@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayFill@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayFill@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySearch@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySearch@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayFind@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayFind@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayMerge@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayMerge@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayPad@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayPad@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayInit@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayInit@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayInsert@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayInsert@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRank@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRank@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRebill@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRebill@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRebind@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRebind@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRefill@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRefill@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRemap@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRemap@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRemove@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRemove@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayReplace@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayReplace@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayResize@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayResize@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayReverse@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayReverse@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySelect@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySelect@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayShuffle@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayShuffle@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySpace@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySpace@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySort@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySort@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySplit@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySplit@routines@split_D1_D0.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arraySplit@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arraySplit@routines@split_D1_D1.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayStrip@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayStrip@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayUnique@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayUnique@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayVerbose@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayVerbose@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_complexAbs@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_complexAbs@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_complexDiv@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_complexDiv@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_complexCompareAll@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAll@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_complexCompareAny@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_complexCompareAny@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_complexCompareLex@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_complexCompareLex@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_cosmicRate@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_cosmicRate@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_cosmology@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_cosmology@routines.inc.F90)
#set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_io@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_io@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distanceEuclid@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distanceEuclid@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distBern@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distBern@routines.inc.F90)


set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distanceMahal@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distanceMahal@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distExp@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distExp@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distGamma@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distGamma@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distGenExpGamma@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distGenExpGamma@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distPiwiPoweto@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distPiwiPoweto@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distPareto@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distPareto@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distPower@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distPower@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distMultiNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distMultiNorm@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines@getUnifRand.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distUnif@routines@setUnifRand.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_distNorm@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_domainBall@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_domainBall@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_except@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_except@routines.inc.F90)


set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_logicalCompare@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_logicalCompare@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_val2logical@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_val2logical@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathCompare@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathCompare@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathCumSum@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathCumSum@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathCumPropExp@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathCumPropExp@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathExp@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathExp@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathFactorial@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathFactorial@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathFactoring@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathFactoring@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathRoot@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathRoot@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathGamma@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathGamma@routines@gamIncLow.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathGamma@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathGamma@routines@gamIncUpp.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathLogAddExp@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathLogAddExp@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathLogSumExp@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathLogSumExp@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_mathMinMax@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_mathMinMax@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_arrayRange@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_arrayRange@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixChol@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixChol@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines@setUpp2Upp.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines@setUppDia2UppDia.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines@setLow2Low.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixCopy@routines@setLowDia2LowDia.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixMulAdd@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixMulAdd@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixMulTri@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixMulTri@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixDet@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixDet@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_matrixInit@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_matrixInit@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_option@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_option@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_val2real@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_val2real@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleCCF@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleCCF@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleCor@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleCor@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleCov@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleCov@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleECDF@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleECDF@routines@ECDF.inc.F90)


set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleMean@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleMean@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleShift@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleShift@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_sampleVar@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_sampleVar@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_timer@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_timer@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/Test_Set_pmod@Subtract_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/Test_Set_pmod@Subtract_smod.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/Test_Set_pmod@SubtractOrdered_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/Test_Set_pmod@SubtractOrdered_smod.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_str@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_str@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_val2str@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_val2str@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_val2int@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_val2int@routines.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_val2complex@routines.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_val2complex@routines.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_partition.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_partition@test_constructPartitionXXXXXX_1.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_partitionRecursiveMaxDen.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_partition.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_partitionRecursiveMinVol.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_partition.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_runSampler_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_runSampler_smod.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_specbase_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_specbase_smod.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_SpecMCMC_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_SpecMCMC_smod.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_specdram_smod.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX@Test_specdram_smod.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraDRAM_ChainFileContents.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_ChainFileContents.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraDISE_ChainFileContents.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_ChainFileContents.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraNest_ChainFileContents.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_ChainFileContents.inc.F90)

set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraDRAM_RefinedChain.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_RefinedChain.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraDISE_RefinedChain.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_RefinedChain.inc.F90)
set_source_files_properties(${paramonte_src_fortran_test_dir}/test_pm_paraNest_RefinedChain.F90 PROPERTIES OBJECT_DEPENDS ${paramonte_src_fortran_test_dir}/test_pm_paraXXXX_RefinedChain.inc.F90)


# Have the .mod files placed in the lib folder
#unset(CMAKE_Fortran_MODULE_DIRECTORY)
#set(CMAKE_Fortran_MODULE_DIRECTORY "${paramonte_bld_test_inc_dir}")

#find_library(PMLIB_PATH ${libname} PATHS ${paramonte_bld_obj_dir})

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Add the test target.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Define the executable name
set(ParaMonteTest_EXE_NAME testParaMonte)

add_executable("${ParaMonteTest_EXE_NAME}" "${ParaMonteTest_SOURCE_FILES}")
set(Test_SRC_BENCHM "${paramonte_src_fortran_test_dir}/Test_PartitionBenchm_pmod.F90")
if ((EXISTS "${SRC_BENCHM}") AND (EXISTS "${Test_SRC_BENCHM}"))
    target_sources("${ParaMonteTest_EXE_NAME}" PRIVATE "${Test_SRC_BENCHM}")
endif()

if (${csid_is_gnu} AND NOT ${codecov_enabled})
    target_link_libraries("${ParaMonteTest_EXE_NAME}" PUBLIC gcov)
endif()

if ("${lib}" MATCHES "[Ss][Tt][Aa][Tt][Ii][Cc]")
    target_link_libraries("${ParaMonteTest_EXE_NAME}" PUBLIC "${libname}")
elseif ("${lib}" MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]" OR "${lib}" MATCHES "[Ss][Hh][Aa][Rr][Ee][Dd]")
    target_link_libraries("${ParaMonteTest_EXE_NAME}" PUBLIC "${libname}")
endif()


set (MACOSX_RPATH TRUE)
#list(APPEND CMAKE_INSTALL_RPATH "${paramonte_bld_lib_dir}")
set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES INSTALL_RPATH "${paramonte_bld_lib_dir}")
set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES LINKER_LANGUAGE "Fortran")
set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${paramonte_bld_test_inc_dir}" )
target_include_directories(${ParaMonteTest_EXE_NAME} PRIVATE "${paramonte_bld_inc_dir}" PRIVATE "${paramonte_bld_test_inc_dir}" )
target_compile_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FC_FLAGS}" "${FCL_FLAGS}" "${fcf}") # "${FPP_FLAGS}"
target_link_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FCL_FLAGS}" "${FL_FLAGS}" "${flf}")

#if ("${test}" STREQUAL "all")
#    target_compile_definitions(${ParaMonteTest_EXE_NAME} PRIVATE "BASIC_TEST_ENABLED=1")
#    target_compile_definitions(${ParaMonteTest_EXE_NAME} PRIVATE "SAMPLER_TEST_ENABLED=1")
#elseif ("${test}" STREQUAL "sampler")
#    target_compile_definitions(${ParaMonteTest_EXE_NAME} PRIVATE "SAMPLER_TEST_ENABLED=1")
#elseif ("${test}" STREQUAL "basic")
#    target_compile_definitions(${ParaMonteTest_EXE_NAME} PRIVATE "BASIC_TEST_ENABLED=1")
#elseif (NOT "${test}" STREQUAL "none")
#    printUsage()
#    message (FATAL_ERROR
#            "\n"
#            "${pmfatal} Unsupported user-specified value for `test`. test=${test}\n"
#            "${pmfatal} Follow the guidelines aboe to properly set the value of `test`.\n"
#            "${pmfatal} Alternatively, skip the testing by dropping the CMake `test` flag.\n"
#            "\n"
#            )
#endif()

# set up BLAS and LAPACK if needed

if (BLAS_FOUND AND LAPACK_FOUND)
    target_link_libraries   (   ${ParaMonteTest_EXE_NAME} PUBLIC
                                ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT}
                            )
endif()

# add parallelization flags if needed

if (OMP_ENABLED)
    set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}" # OpenMP_Fortran_FLAGS is a CMake variable.
                          LINK_FLAGS    "${OpenMP_Fortran_FLAGS}"
                          )
endif()

if (MPI_ENABLED)
    set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    target_include_directories(${ParaMonteTest_EXE_NAME} PUBLIC "${MPI_Fortran_INCLUDE_PATH}")
    target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC "${MPI_Fortran_LIBRARIES}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Report the test executable build options.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

get_target_property(PMLIB_TEST_INCLUDE_DIRECTORIES   ${ParaMonteTest_EXE_NAME} INCLUDE_DIRECTORIES)
get_target_property(PMLIB_TEST_COMPILE_DEFINITIONS   ${ParaMonteTest_EXE_NAME} COMPILE_DEFINITIONS)
get_target_property(PMLIB_TEST_COMPILE_OPTIONS       ${ParaMonteTest_EXE_NAME} COMPILE_OPTIONS)
get_target_property(PMLIB_TEST_COMPILE_FLAGS         ${ParaMonteTest_EXE_NAME} COMPILE_FLAGS)

get_target_property(PMLIB_TEST_LINK_DIRECTORIES  ${ParaMonteTest_EXE_NAME} LINK_DIRECTORIES)
get_target_property(PMLIB_TEST_INSTALL_NAME_DIR  ${ParaMonteTest_EXE_NAME} INSTALL_NAME_DIR)
get_target_property(PMLIB_TEST_LINK_OPTIONS      ${ParaMonteTest_EXE_NAME} LINK_OPTIONS)
get_target_property(PMLIB_TEST_LINK_FLAGS        ${ParaMonteTest_EXE_NAME} LINK_FLAGS)

message(NOTICE "${pmattn} ParaMonte library test target include directories : ${PMLIB_TEST_INCLUDE_DIRECTORIES}" )
message(NOTICE "${pmattn} ParaMonte library test target compile definitions : ${PMLIB_TEST_COMPILE_DEFINITIONS}" )
message(NOTICE "${pmattn} ParaMonte library test target compile options     : ${PMLIB_TEST_COMPILE_OPTIONS}" )
message(NOTICE "${pmattn} ParaMonte library test target compile flags       : ${PMLIB_TEST_COMPILE_FLAGS}" )
message(NOTICE "${pmattn} ParaMonte library test target link directories    : ${PMLIB_TEST_LINK_DIRECTORIES}" )
message(NOTICE "${pmattn} ParaMonte library test target install dir name    : ${PMLIB_TEST_INSTALL_NAME_DIR}" )
message(NOTICE "${pmattn} ParaMonte library test target link options        : ${PMLIB_TEST_LINK_OPTIONS}" )
message(NOTICE "${pmattn} ParaMonte library test target link flags          : ${PMLIB_TEST_LINK_FLAGS}" )
message(NOTICE "${pmattn} ParaMonte library test rapth                      : ${CMAKE_INSTALL_RPATH}" )
message(NOTICE "${pmattn} paramonte_bld_lib_dir                             : ${paramonte_bld_lib_dir}" )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Install the test executable.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(CMAKE_INSTALL_PREFIX "${paramonte_bld_test_bin_dir}/")
install(TARGETS ${ParaMonteTest_EXE_NAME}
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}"
        )

#if (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]")
#    set(CMAKE_INSTALL_PREFIX "${paramonte_bld_test_bin_dir}/")
#    install(TARGETS ${libname}
#            CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
#            DESTINATION "${CMAKE_INSTALL_PREFIX}"
#            )
#endif()

message(NOTICE "${pmattn} CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )

#if (NOT CMAKE_FORTRAN_LINK_EXECUTABLE)
#    set(CMAKE_FORTRAN_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_FLAGS> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
#    message(NOTICE "CMAKE_FORTRAN_LINK_EXECUTABLE: ${CMAKE_FORTRAN_LINK_EXECUTABLE}")
#endif ()

# copy input data for testing to the _bin folder

# enable_testing() # only useful when ctests is used, which is not the case here.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Add rule for copying the test input files.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

message(NOTICE "${pmattn} Copying the tests input files" )
message(NOTICE "${pmattn} from: ${paramonte_src_fortran_test_input_dir}/" )
message(NOTICE "${pmattn}   to: ${paramonte_bld_test_bin_input_dir}/" )

add_custom_target(  copyTestInputFiles ALL
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    "${paramonte_src_fortran_test_input_dir}"
                    "${paramonte_bld_test_bin_input_dir}"
                    )
add_dependencies(${ParaMonteTest_EXE_NAME} copyTestInputFiles)

## ${TEST_INPUT_FILES} is a semicolon-separated list and should not be quoted below.

#file(GLOB TEST_INPUT_FILES "${paramonte_src_fortran_test_input_dir}/*")
#foreach(input_file ${TEST_INPUT_FILES})
#    add_custom_command( TARGET ${ParaMonteTest_EXE_NAME} POST_BUILD
#                        # ${CMAKE_COMMAND} -E achieves a cross-platform copy
#                        COMMAND ${CMAKE_COMMAND} -E copy
#                        "${input_file}"
#                        "${paramonte_bld_test_bin_input_dir}/"
#                        COMMENT "${pmattn} Copying test data ${input_file} to the test runtime directory ..."
#                        # $<TARGET_FILE:${ParaMonteTest_EXE_NAME}> resolves to the full path of the library created by ParaMonteTest_EXE_NAME: "${paramonte_bld_test_bin_dir}/", regardless of platform or compiler
#                        )
#endforeach(input_file)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: Add rule for running the tests.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

add_custom_target(  test
                    COMMAND "${paramonte_bld_test_bin_dir}/${ParaMonteTest_EXE_NAME}"
                    )
#add_dependencies(${ParaMonteTest_EXE_NAME} runParaMonteTest)
