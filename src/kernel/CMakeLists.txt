####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a
####   copy of this software and associated documentation files (the "Software"),
####   to deal in the Software without restriction, including without limitation
####   the rights to use, copy, modify, merge, publish, distribute, sublicense,
####   and/or sell copies of the Software, and to permit persons to whom the
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your
####   work (education/research/industry/development/...) by citing the ParaMonte
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

# ParaMonte library
SET(ParaMonte_SOURCE_FILES
   #${ParaMonteKernel_SRC_DIR}/ifport.f90
    ${ParaMonteKernel_SRC_DIR}/Constants_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Misc_mod.f90
   #${ParaMonteKernel_SRC_DIR}/EconomicsToolbox_mod.f90
    ${ParaMonteKernel_SRC_DIR}/JaggedArray_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Decoration_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Decoration_mod@Routines_smod.f90
    ${ParaMonteKernel_SRC_DIR}/Err_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Optimization_mod.f90
    ${ParaMonteKernel_SRC_DIR}/DateTime_mod.f90
    ${ParaMonteKernel_SRC_DIR}/String_mod.f90
    ${ParaMonteKernel_SRC_DIR}/FileContents_mod.f90
    ${ParaMonteKernel_SRC_DIR}/System_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Path_mod.f90
    ${ParaMonteKernel_SRC_DIR}/FileList_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Timer_mod.f90
    ${ParaMonteKernel_SRC_DIR}/TimerCPU_mod.f90
    ${ParaMonteKernel_SRC_DIR}/RandomSeed_mod.f90
    ${ParaMonteKernel_SRC_DIR}/File_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Sort_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Math_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Matrix_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Statistics_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Kmeans_mod.f90
    ${ParaMonteKernel_SRC_DIR}/KmeansOOP_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Partition_mod.f90
    ${ParaMonteKernel_SRC_DIR}/PartitionMaxDen_mod.f90
    ${ParaMonteKernel_SRC_DIR}/PartitionMinVol_mod.f90
    ${ParaMonteKernel_SRC_DIR}/GeoCyclicFit_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Parallelism_mod.f90
    ${ParaMonteKernel_SRC_DIR}/CorrCoef_mod.f90
    ${ParaMonteKernel_SRC_DIR}/TranGaus_mod.f90
    ${ParaMonteKernel_SRC_DIR}/CrossCorr_mod.f90
    ${ParaMonteKernel_SRC_DIR}/Integration_mod.f90
    ${ParaMonteKernel_SRC_DIR}/QuadPackSPR_mod.f90
    #${ParaMonteKernel_SRC_DIR}/QuadPackDPR_mod.f
    ${ParaMonteKernel_SRC_DIR}/Cosmology_mod.f90
    ${ParaMonteKernel_SRC_DIR}/StarFormation_mod.f90
    ${ParaMonteKernel_SRC_DIR}/BandSpectrum_mod.f90
    ${ParaMonteKernel_SRC_DIR}/BandSpectrum_mod@PhotonFluence_smod.f90
    ${ParaMonteKernel_SRC_DIR}/BandSpectrum_mod@EnergyFluence_smod.f90
    ${ParaMonteKernel_SRC_DIR}/Batse_mod.f90
    # SpecBase
    ${ParaMonteKernel_SRC_DIR}/SpecBase_SampleSize_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_RandomSeed_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_Description_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_InterfaceType_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_OutputFileName_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_OutputDelimiter_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_ChainFileFormat_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_VariableNameList_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_RestartFileFormat_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_OutputColumnWidth_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_OverwriteRequested_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_SystemInfoFilePath_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_OutputRealPrecision_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_SilentModeRequested_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_DomainLowerLimitVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_DomainUpperLimitVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_ParallelizationModel_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_InputFileHasPriority_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_ProgressReportPeriod_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_TargetAcceptanceRate_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_MpiFinalizeRequested_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecBase_mod.f90
    # SpecMCMC
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ChainSize_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ScaleFactor_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_StartPointVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalModel_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCovMat_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCorMat_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartStdVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementCount_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementMethod_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointRequested_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecMCMC_mod.f90
    # SpecDRAM
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionCount_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecDRAM_mod.f90
    # SpecNest
    ${ParaMonteKernel_SRC_DIR}/SpecNest_Tightness_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_Tolerance_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_ScaleFactor_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_ProposalModel_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_LiveSampleSize_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_InclusionFraction_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_AdaptiveUpdateCount_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_AdaptiveUpdatePeriod_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_MahalSqWeightExponent_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_StabilizationRequested_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_MaxAllowedKmeansFailure_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_MaxKvolumeLoopRecursion_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_MaxAllowedMinVolFailure_mod.f90
    ${ParaMonteKernel_SRC_DIR}/SpecNest_mod.f90
    # Samplers
    ${ParaMonteKernel_SRC_DIR}/ParaMonte_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaMonteLogFunc_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ChainFileContents_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_ChainFileContents_mod.f90
    #${ParaMonteKernel_SRC_DIR}/ParaNest_ChainFileContents_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_RefinedChain_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_RefinedChain_mod.f90
    #${ParaMonteKernel_SRC_DIR}/ParaNest_RefinedChain_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalAbstract_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalAbstract_mod.f90
    #${ParaMonteKernel_SRC_DIR}/ParaNest_ProposalAbstract_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalUniform_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalUniform_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalNormal_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalNormal_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaMCMC_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_mod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Setup_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Setup_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Input_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Input_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Kernel_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Postprocess_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Postprocess_smod.f90
    ${ParaMonteKernel_SRC_DIR}/ParaMonte@runParaDRAM.f90
)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/PartitionMaxDen_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/Partition_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/PartitionMinVol_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/Partition_mod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaMonte_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaMonte_mod@version@kernel.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_ChainFileContents_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_ChainFileContents_mod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Input_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Input_smod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Input_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Input_smod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Setup_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Setup_smod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Setup_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Setup_smod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Kernel_smod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Kernel_smod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Kernel_smod@nextMove.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Kernel_smod@nextMove.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Postprocess_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Postprocess_smod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Postprocess_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod@Postprocess_smod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalNormal_mod.f90    PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_Proposal_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalNormal_mod.f90    PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_Proposal_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalNormal_mod.f90    PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalAbstract_mod.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalNormal_mod.f90    PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalAbstract_mod.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalUniform_mod.f90   PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_Proposal_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalUniform_mod.f90   PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_Proposal_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalUniform_mod.f90   PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalAbstract_mod.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalUniform_mod.f90   PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalAbstract_mod.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ProposalAbstract_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_ProposalAbstract_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ProposalAbstract_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_ProposalAbstract_mod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_RefinedChain_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_RefinedChain_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_RefinedChain_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_RefinedChain_mod.inc.f90)

set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_ChainFileContents_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod.inc.f90)
set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_ChainFileContents_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteKernel_SRC_DIR}/ParaXXXX_mod.inc.f90)

#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_SampleSize_mod.f90                          PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_SampleSize_desc.f90                          ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_RandomSeed_mod.f90                          PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_RandomSeed_desc.f90                          ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_Description_mod.f90                         PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_Description_desc.f90                         ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_InterfaceType_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_InterfaceType_desc.f90                       ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_OutputFileName_mod.f90                      PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_OutputFileName_desc.f90                      ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_OutputDelimiter_mod.f90                     PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_OutputDelimiter_desc.f90                     ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_ChainFileFormat_mod.f90                     PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_ChainFileFormat_desc.f90                     ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_VariableNameList_mod.f90                    PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_VariableNameList_desc.f90                    ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_RestartFileFormat_mod.f90                   PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_RestartFileFormat_desc.f90                   ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_OutputColumnWidth_mod.f90                   PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_OutputColumnWidth_desc.f90                   ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_OutputRealPrecision_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_OutputRealPrecision_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_SilentModeRequested_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_SilentModeRequested_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_DomainLowerLimitVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_DomainLowerLimitVec_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_DomainUpperLimitVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_DomainUpperLimitVec_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_ParallelizationModel_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_ParallelizationModel_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_InputFileHasPriority_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_InputFileHasPriority_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_ProgressReportPeriod_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_ProgressReportPeriod_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_TargetAcceptanceRate_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_TargetAcceptanceRate_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_MpiFinalizeRequested_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_MpiFinalizeRequested_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_desc.f90             ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_desc.f90             ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ChainSize_mod.f90                           PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ChainSize_desc.f90                           ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ScaleFactor_mod.f90                         PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ScaleFactor_desc.f90                         ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_StartPointVec_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_StartPointVec_desc.f90                       ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalModel_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalModel_desc.f90                       ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCovMat_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCovMat_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCorMat_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartCorMat_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartStdVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_ProposalStartStdVec_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementCount_desc.f90               ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementMethod_mod.f90              PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_SampleRefinementMethod_desc.f90              ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointRequested_mod.f90           PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointRequested_desc.f90           ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_mod.f90 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_desc.f90 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_mod.f90 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_desc.f90 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_desc.f90                 ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_desc.f90                ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_desc.f90               ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionCount_desc.f90               ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_desc.f90             ")
#set_source_files_properties(${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_mod.f90      PROPERTIES OBJECT_DEPENDS "${ParaMonteKernel_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_desc.f90      ")

# enable debug info on IO source files

#if (gnu_compiler)
#    set(SRC_LEVEL_DEBUG_FLAG "-g")
#elseif(intel_compiler)
#    if (WIN32)
#        set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Input_smod.f90 PROPERTIES COMPILE_FLAGS "/debug")
#        set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Input_smod.f90 PROPERTIES COMPILE_FLAGS "/debug")
#    elseif (UNIX OR APPLE)
#        set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDRAM_mod@Input_smod.f90 PROPERTIES COMPILE_FLAGS "-g")
#        set_source_files_properties(${ParaMonteKernel_SRC_DIR}/ParaDISE_mod@Input_smod.f90 PROPERTIES COMPILE_FLAGS "-g") # ${INTEL_Fortran_DEBUG_FLAGS}
#    endif()
#endif()

# search for static libraries

function(find_static_library LIB_NAME OUT)

    if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    elseif (UNIX)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()

    find_library(
        FOUND_${LIB_NAME}_STATIC
        ${LIB_NAME}
        )

    if (FOUND_${LIB_NAME}_STATIC)
        get_filename_component(ABS_FILE ${FOUND_${LIB_NAME}_STATIC} ABSOLUTE)
        set(${OUT} ${ABS_FILE} PARENT_SCOPE)
    else()
        message ( WARNING
                " \n"
                " ${pmwarn}\n"
                " ${pmattn} Unable to find ${LIB_NAME}${CMAKE_FIND_LIBRARY_SUFFIXES}\n"
                " ${pmattn} Please report any failures to shahmoradi@utexas.edu"
                )
        set(${OUT} "" PARENT_SCOPE)
    endif()

endfunction()

unset(FL_FLAGS_PMLIB)
unset(FCL_FLAGS_PMLIB)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# set shared library Fortran linker flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (pmlib_shared)

    # The following does not seem to work.

    #if (gnu_compiler AND UNIX AND NOT APPLE)
    #    #set(CMAKE_SKIP_RPATH ON CACHE BOOL "Skip adding a build-dir rpath")
    #    set(CMAKE_BUILD_WITH_INSTALL_RPATH ON CACHE BOOL "Build using the install rpath")
    #    set(CMAKE_BUILD_RPATH_USE_ORIGIN ON CACHE BOOL "Use relative rpaths")
    #    set(CMAKE_INSTALL_RPATH "$ORIGIN")
    #endif()

    add_library(${PMLIB_NAME} SHARED ${ParaMonte_SOURCE_FILES})

    #if (gnu_compiler)
        #set_target_properties( ${PMLIB_NAME} PROPERTIES COMPILE_FLAGS "-static-libgfortran" )
        #set_target_properties( ${PMLIB_NAME} PROPERTIES LINK_FLAGS "-static-libgfortran" )
        #set_target_properties( ${PMLIB_NAME} PROPERTIES COMPILE_FLAGS "-Wl,--whole-archive -lgfortran -static -lgfortran" )
        #set(CMAKE_SHARED_LINKER_FLAGS "-static-libgfortran")
        #set_target_properties( ${PMLIB_NAME} PROPERTIES LINK_FLAGS "-static-libgfortran" )
        #if (APPLE)
            #set_target_properties( ${PMLIB_NAME} PROPERTIES LINK_FLAGS "-Wl,--whole-archive" )
        #elseif(UNIX)
            #set_target_properties( ${PMLIB_NAME} PROPERTIES LINK_FLAGS "-Wl,-static-libgfortran" )
        #endif()

        #set_target_properties( ${PMLIB_NAME} PROPERTIES SKIP_BUILD_RPATH FALSE) # don't skip the full RPATH for the build tree
        #set_target_properties( ${PMLIB_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE) # when building, don't use the install RPATH for the build
        #set_target_properties( ${PMLIB_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE) # if TRUE, append the automatically determined parts of RPATH
        #set_target_properties( ${PMLIB_NAME} PROPERTIES INSTALL_RPATH "./") # the RPATH to be used when installing the library
        #set_target_properties( ${PMLIB_NAME} PROPERTIES LINK_FLAGS "-Wl,-rpath,'$ORIGIN'" ) # if uncommented, this will override all of the above!
    #endif()

   #if (gnu_compiler)
   #    find_static_library(libgfortran LIBGFORTRAN_PATH)
   #    find_static_library(libquadmath LIBQUADMATH_PATH)
   #    if (${BTYPE} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries(${PMLIB_NAME} debug ${LIBGFORTRAN_PATH})
   #        target_link_libraries(${PMLIB_NAME} debug ${LIBQUADMATH_PATH})
   #    else() #if (${BTYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries(${PMLIB_NAME} optimized ${LIBGFORTRAN_PATH})
   #        target_link_libraries(${PMLIB_NAME} optimized ${LIBQUADMATH_PATH})
   #    endif()
   #elseif(intel_compiler)
   #    find_static_library(libifport LIBIFPORT_PATH)
   #    if (${BTYPE} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries(${PMLIB_NAME} debug ${LIBIFPORT_PATH})
   #    else() #if (${BTYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries(${PMLIB_NAME} optimized ${LIBIFPORT_PATH})
   #    endif()
   #endif()

    set_target_properties(${PMLIB_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
    set_target_properties(${PMLIB_NAME} PROPERTIES BUILD_RPATH_USE_ORIGIN ON)

    unset(RPATH_PROP)
    if (APPLE)
        set(RPATH_PROP "@loader_path")
    elseif(UNIX)
        set(RPATH_PROP "$ORIGIN")
    endif()

    if (DEFINED MPILIB_NAME)
        set(RPATH_PROP "${RPATH_PROP}:${RPATH_PROP}/${MPILIB_NAME}")
    endif()

    if (DEFINED RPATH_PROP)
        set_target_properties(${PMLIB_NAME} PROPERTIES INSTALL_RPATH "${RPATH_PROP}")
    endif()

    #### set up the library-specific compile link flags

    if(gnu_compiler)

        set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
        -shared
        -fPIC
        )

        set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
        -shared
        )

        if (deploy_enabled)
        #    set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
        #    -static-libgfortran
        #    -static-libgcc
        #    )
        endif()

    elseif(intel_compiler)

        if (WIN32)

            set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
            #/threads " # these flags are actually included by default in recent ifort implementations
            /libs:dll
            )

            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            /dll
            )

        elseif(UNIX)

            set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
            -fpic # Request compiler to generate position-independent code.
            )

            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            -qopenmp-link=static
            -static-intel
            )

            if (deploy_enabled)
                set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
                -static-libstdc++
                -static-libgcc
                )
            endif()

            if(APPLE)

                set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
                -noall_load
                # -weak_references_mismatches non-weak -threads -arch_only i386
                )

                set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
                -noall_load
                -dynamiclib
                )

            else()

                set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
                -shared
                )

            endif()

        endif()

        #### Include MKL

        #if (CMAKE_BUILD_TYPE MATCHES "Release|RELEASE|release")
        #    #find_package(MKL REQUIRED)
        #    target_link_directories(${PMLIB_NAME} PUBLIC ${MKL_LIBRARIES})
        #    target_include_directories(${PMLIB_NAME} PUBLIC ${MKL_INCLUDE_DIRS})
        #    target_link_libraries(${PMLIB_NAME} PUBLIC mkl_intel_lp64 mkl_sequential mkl_core)
        #endif()

    endif()

else() # static lib

    add_library(${PMLIB_NAME} STATIC ${ParaMonte_SOURCE_FILES})

endif()

if (MPI_ENABLED)
    #include_directories(${MPI_Fortran_INCLUDE_PATH})
    target_include_directories(${PMLIB_NAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(${PMLIB_NAME} PUBLIC ${MPI_Fortran_LIBRARIES})
endif()

#if (NOT WIN32)
#    if (DEFINED MATLAB_ROOT_DIR)
#        target_include_directories(${PMLIB_NAME} PUBLIC "${MATLAB_ROOT_DIR}/extern/include")
#        target_link_directories(${PMLIB_NAME} BEFORE
#        PUBLIC "$(MATLAB_ROOT_DIR)/bin/${MATLAB_PLATFORM_SUBDIR}"
#        PUBLIC "$(MATLAB_ROOT_DIR)/extern/bin/${MATLAB_PLATFORM_SUBDIR}"
#        )
#        target_link_libraries(${PMLIB_NAME} PUBLIC mx PUBLIC mex PUBLIC mat PUBLIC m PUBLIC stdc++)
#    endif()
#endif()

# ensure preprocessor flags and module dependencies containing preprocessor flags are respected

# define the output modules directory

set_target_properties( ${PMLIB_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${ParaMonte_MOD_DIR}" )
# target_compile_definitions( ${PMLIB_NAME} INTERFACE "${FPP_FLAGS}" ) # cmake currently does not add preprocessor automatically

#string(REPLACE ";" " " FCL_FLAGS ${FCL_FLAGS})
#message(STATUS "${pmattn} ParaMonte library build flags: ${FPP_FLAGS} ${FCL_FLAGS}")
target_compile_options(${PMLIB_NAME} BEFORE PRIVATE "${FPP_FLAGS}" "${FC_FLAGS}" "${FCL_FLAGS}" "${FCL_FLAGS_PMLIB}")

#message( STATUS "${pmattn} ParaMonte library build linker flags: ${${FCL_FLAGS}${FL_FLAGS_PMLIB}}${MATLAB_LINKER_FLAGS}" )
#target_link_options( ${PMLIB_NAME} BEFORE PRIVATE ${${FCL_FLAGS}${FL_FLAGS_PMLIB}}${MATLAB_LINKER_FLAGS} )
#message( STATUS "${pmattn} ParaMonte library build linker flags: ${FCL_FLAGS} ${FL_FLAGS_PMLIB}" )
#target_link_options( ${PMLIB_NAME} BEFORE PRIVATE ${${FCL_FLAGS}${FL_FLAGS}${FL_FLAGS_PMLIB}} )
target_link_options( ${PMLIB_NAME} BEFORE PRIVATE "${FCL_FLAGS}" "${FL_FLAGS}" "${FL_FLAGS_PMLIB}" )

get_target_property(PMLIB_INCLUDE_DIRECTORIES   ${PMLIB_NAME} INCLUDE_DIRECTORIES)
get_target_property(PMLIB_COMPILE_DEFINITIONS   ${PMLIB_NAME} COMPILE_DEFINITIONS)
get_target_property(PMLIB_COMPILE_OPTIONS       ${PMLIB_NAME} COMPILE_OPTIONS)
get_target_property(PMLIB_COMPILE_FLAGS         ${PMLIB_NAME} COMPILE_FLAGS)

get_target_property(PMLIB_LINK_DIRECTORIES  ${PMLIB_NAME} LINK_DIRECTORIES)
get_target_property(PMLIB_INSTALL_NAME_DIR  ${PMLIB_NAME} INSTALL_NAME_DIR)
get_target_property(PMLIB_LINK_OPTIONS      ${PMLIB_NAME} LINK_OPTIONS)
get_target_property(PMLIB_LINK_FLAGS        ${PMLIB_NAME} LINK_FLAGS)

message( STATUS "${pmattn} target include directories: ${PMLIB_INCLUDE_DIRECTORIES}" )
message( STATUS "${pmattn} target compile definitions: ${PMLIB_COMPILE_DEFINITIONS}" )
message( STATUS "${pmattn} target compile options    : ${PMLIB_COMPILE_OPTIONS}" )
message( STATUS "${pmattn} target compile flags      : ${PMLIB_COMPILE_FLAGS}" )
message( STATUS "${pmattn} target link directories   : ${PMLIB_LINK_DIRECTORIES}" )
message( STATUS "${pmattn} target install name dir   : ${PMLIB_INSTALL_NAME_DIR}" )
message( STATUS "${pmattn} target link options       : ${PMLIB_LINK_OPTIONS}" )
message( STATUS "${pmattn} target link flags         : ${PMLIB_LINK_FLAGS}" )
message( STATUS "${pmattn} CMAKE INSTALL PATH        : ${ParaMonte_LIB_DIR}" )

if (CODECOV_ENABLED AND gnu_compiler)
    target_link_libraries(${PMLIB_NAME} PUBLIC gcov)
endif()

install(TARGETS ${PMLIB_NAME}
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        DESTINATION "${ParaMonte_LIB_DIR}/"
        )
