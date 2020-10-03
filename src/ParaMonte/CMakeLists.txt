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
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

# ParaMonte library
SET(ParaMonte_SOURCE_FILES
   #${ParaMonte_SRC_DIR}/ifport.f90
    ${ParaMonte_SRC_DIR}/Constants_mod.f90
    ${ParaMonte_SRC_DIR}/Misc_mod.f90
   #${ParaMonte_SRC_DIR}/EconomicsToolbox_mod.f90
    ${ParaMonte_SRC_DIR}/JaggedArray_mod.f90
    ${ParaMonte_SRC_DIR}/Decoration_mod.f90
    ${ParaMonte_SRC_DIR}/Decoration_mod@Routines_smod.f90
    ${ParaMonte_SRC_DIR}/Err_mod.f90
    ${ParaMonte_SRC_DIR}/Optimization_mod.f90
    ${ParaMonte_SRC_DIR}/DateTime_mod.f90
    ${ParaMonte_SRC_DIR}/String_mod.f90
    ${ParaMonte_SRC_DIR}/System_mod.f90
    ${ParaMonte_SRC_DIR}/Path_mod.f90
    ${ParaMonte_SRC_DIR}/FileList_mod.f90
    ${ParaMonte_SRC_DIR}/Timer_mod.f90
    ${ParaMonte_SRC_DIR}/TimerCPU_mod.f90
    ${ParaMonte_SRC_DIR}/RandomSeed_mod.f90
    ${ParaMonte_SRC_DIR}/FileContents_mod.f90
    ${ParaMonte_SRC_DIR}/File_mod.f90
    ${ParaMonte_SRC_DIR}/Sort_mod.f90
    ${ParaMonte_SRC_DIR}/Math_mod.f90
    ${ParaMonte_SRC_DIR}/Matrix_mod.f90
    ${ParaMonte_SRC_DIR}/Statistics_mod.f90
    ${ParaMonte_SRC_DIR}/Parallelism_mod.f90
    ${ParaMonte_SRC_DIR}/CorrCoef_mod.f90
    ${ParaMonte_SRC_DIR}/TranGaus_mod.f90
    ${ParaMonte_SRC_DIR}/CrossCorr_mod.f90
    ${ParaMonte_SRC_DIR}/Integration_mod.f90
    ${ParaMonte_SRC_DIR}/QuadPackSPR_mod.f90
    ${ParaMonte_SRC_DIR}/QuadPackDPR_mod.f
    ${ParaMonte_SRC_DIR}/Cosmology_mod.f90
    ${ParaMonte_SRC_DIR}/StarFormation_mod.f90
    ${ParaMonte_SRC_DIR}/BandSpectrum_mod.f90
    ${ParaMonte_SRC_DIR}/Batse_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_SampleSize_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_RandomSeed_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_Description_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_InterfaceType_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_OutputFileName_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_OutputDelimiter_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_ChainFileFormat_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_VariableNameList_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_RestartFileFormat_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_OutputColumnWidth_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_OverwriteRequested_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_SystemInfoFilePath_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_OutputRealPrecision_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_SilentModeRequested_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_DomainLowerLimitVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_DomainUpperLimitVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_ParallelizationModel_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_InputFileHasPriority_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_ProgressReportPeriod_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_TargetAcceptanceRate_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_MpiFinalizeRequested_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_mod.f90
    ${ParaMonte_SRC_DIR}/SpecBase_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ChainSize_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ScaleFactor_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_StartPointVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ProposalModel_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCovMat_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCorMat_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartStdVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementCount_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementMethod_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointRequested_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecMCMC_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionCount_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_mod.f90
    ${ParaMonte_SRC_DIR}/SpecDRAM_mod.f90
    ${ParaMonte_SRC_DIR}/ParaMonteLogFunc_mod.f90
    ${ParaMonte_SRC_DIR}/ParaMonteChainFileContents_mod.f90
    ${ParaMonte_SRC_DIR}/ParaMCMCRefinedChain_mod.f90
    ${ParaMonte_SRC_DIR}/ParaMonte_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDISEProposalAbstract_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAMProposalAbstract_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAMProposalNormal_mod.f90
    ${ParaMonte_SRC_DIR}/ParaMCMC_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAM_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDISE_mod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90
    ${ParaMonte_SRC_DIR}/ParaDISE_mod@Kernel_smod.f90
    ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Setup_smod.f90
    ${ParaMonte_SRC_DIR}/ParaDISE_mod@Setup_smod.f90
    ${ParaMonte_SRC_DIR}/ParaMonte@runParaDRAM.f90
)
   #${ParaMonte_SRC_DIR}/ParaDRAMProposalSymmetricUniform_mod.f90

set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAM_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAM_mod@Setup_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Setup_smod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDISE_mod@Setup_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Setup_smod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod@nextMove.inc.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAM_mod@Kernel_smod@markovGen.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalNormal_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAMProposal_mod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAMProposal_mod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalAbstract_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAMProposalAbstract_mod.inc.f90)
set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDISEProposalAbstract_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDRAMProposalAbstract_mod.inc.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaDISEProposalAbstract_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaMonteLogFunc_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/Statistics_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/ParaMonte_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/Constants_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/SpecBase_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/SpecMCMC_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/SpecDRAM_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/String_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/Matrix_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/Math_mod.f90)
#set_source_files_properties(${ParaMonte_SRC_DIR}/ParaDRAMProposalUniform_mod.f90  PROPERTIES OBJECT_DEPENDS ${ParaMonte_SRC_DIR}/Err_mod.f90)

#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_SampleSize_mod.f90                          PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_SampleSize_desc.f90                          ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_RandomSeed_mod.f90                          PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_RandomSeed_desc.f90                          ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_Description_mod.f90                         PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_Description_desc.f90                         ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_InterfaceType_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_InterfaceType_desc.f90                       ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_OutputFileName_mod.f90                      PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_OutputFileName_desc.f90                      ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_OutputDelimiter_mod.f90                     PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_OutputDelimiter_desc.f90                     ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_ChainFileFormat_mod.f90                     PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_ChainFileFormat_desc.f90                     ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_VariableNameList_mod.f90                    PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_VariableNameList_desc.f90                    ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_RestartFileFormat_mod.f90                   PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_RestartFileFormat_desc.f90                   ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_OutputColumnWidth_mod.f90                   PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_OutputColumnWidth_desc.f90                   ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_OutputRealPrecision_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_OutputRealPrecision_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_SilentModeRequested_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_SilentModeRequested_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_DomainLowerLimitVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_DomainLowerLimitVec_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_DomainUpperLimitVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_DomainUpperLimitVec_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_ParallelizationModel_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_ParallelizationModel_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_InputFileHasPriority_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_InputFileHasPriority_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_ProgressReportPeriod_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_ProgressReportPeriod_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_TargetAcceptanceRate_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_TargetAcceptanceRate_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_MpiFinalizeRequested_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_MpiFinalizeRequested_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToWarn_desc.f90             ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecBase_MaxNumDomainCheckToStop_desc.f90             ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ChainSize_mod.f90                           PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ChainSize_desc.f90                           ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ScaleFactor_mod.f90                         PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ScaleFactor_desc.f90                         ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_StartPointVec_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_StartPointVec_desc.f90                       ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ProposalModel_mod.f90                       PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ProposalModel_desc.f90                       ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCovMat_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCovMat_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCorMat_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartCorMat_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartStdVec_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_ProposalStartStdVec_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementCount_desc.f90               ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementMethod_mod.f90              PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_SampleRefinementMethod_desc.f90              ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointRequested_mod.f90           PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointRequested_desc.f90           ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_mod.f90 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainLowerLimitVec_desc.f90 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_mod.f90 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecMCMC_RandomStartPointDomainUpperLimitVec_desc.f90 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_mod.f90                 PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdateCount_desc.f90                 ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_mod.f90                PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_AdaptiveUpdatePeriod_desc.f90                ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_GreedyAdaptationCount_desc.f90               ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionCount_mod.f90               PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionCount_desc.f90               ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_mod.f90             PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_BurninAdaptationMeasure_desc.f90             ")
#set_source_files_properties(${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_mod.f90      PROPERTIES OBJECT_DEPENDS "${ParaMonte_SRC_DIR}/SpecDRAM_DelayedRejectionScaleFactorVec_desc.f90      ")

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

# Have the .mod files placed in the lib folder

if (${LTYPE} MATCHES "[Ss][Tt][Aa][Tt][Ii][Cc]")

    add_library(${PMLIB} STATIC ${ParaMonte_SOURCE_FILES})

elseif (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]")

    add_library(${PMLIB} SHARED ${ParaMonte_SOURCE_FILES})
   #add_library(${PMLIB} MODULE ${ParaMonte_SOURCE_FILES})

    if (gnu_compiler)
        #set_target_properties( ${PMLIB} PROPERTIES COMPILE_FLAGS "-static-libgfortran" )
        #set_target_properties( ${PMLIB} PROPERTIES LINK_FLAGS "-static-libgfortran" )
        #set_target_properties( ${PMLIB} PROPERTIES COMPILE_FLAGS "-Wl,--whole-archive -lgfortran -static -lgfortran" )
        #set(CMAKE_SHARED_LINKER_FLAGS "-static-libgfortran")
        #set_target_properties( ${PMLIB} PROPERTIES LINK_FLAGS "-static-libgfortran" )
        #if (APPLE)
            #set_target_properties( ${PMLIB} PROPERTIES LINK_FLAGS "-Wl,--whole-archive" )
        #elseif(UNIX)
            #set_target_properties( ${PMLIB} PROPERTIES LINK_FLAGS "-Wl,-static-libgfortran" )
        #endif()
    elseif(intel_compiler)
        if (UNIX OR APPLE)
            set_target_properties( ${PMLIB} PROPERTIES LINK_FLAGS "-static-intel" )
        endif()
    endif()

   #if (gnu_compiler)
   #    find_static_library(libgfortran LIBGFORTRAN_PATH)
   #    find_static_library(libquadmath LIBQUADMATH_PATH)
   #    if (${BTYPE} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries(${PMLIB} debug ${LIBGFORTRAN_PATH})
   #        target_link_libraries(${PMLIB} debug ${LIBQUADMATH_PATH})
   #    else() #if (${BTYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries(${PMLIB} optimized ${LIBGFORTRAN_PATH})
   #        target_link_libraries(${PMLIB} optimized ${LIBQUADMATH_PATH})
   #    endif()
   #elseif(intel_compiler)
   #    find_static_library(libifport LIBIFPORT_PATH)
   #    if (${BTYPE} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries(${PMLIB} debug ${LIBIFPORT_PATH})
   #    else() #if (${BTYPE} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries(${PMLIB} optimized ${LIBIFPORT_PATH})
   #    endif()
   #endif()

else()

    message ( FATAL_ERROR 
            " \n"
            "${pmfatal}\n"
            "   Unrecognized library type specified for building ParaMonte library:\n"
            " \n"
            "   LTYPE=${LTYPE}\n"
            " \n"
            )    

endif()

if (MPI_ENABLED)
    #include_directories(${MPI_Fortran_INCLUDE_PATH})
    target_include_directories(${PMLIB} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(${PMLIB} PUBLIC ${MPI_Fortran_LIBRARIES})
endif()

#if (NOT WIN32)
#    if (DEFINED MATLAB_ROOT_DIR)
#        target_include_directories(${PMLIB} PUBLIC "${MATLAB_ROOT_DIR}/extern/include")
#        target_link_directories(${PMLIB} BEFORE 
#        PUBLIC "$(MATLAB_ROOT_DIR)/bin/${MATLAB_PLATFORM_SUBDIR}"
#        PUBLIC "$(MATLAB_ROOT_DIR)/extern/bin/${MATLAB_PLATFORM_SUBDIR}"
#        )
#        target_link_libraries(${PMLIB} PUBLIC mx PUBLIC mex PUBLIC mat PUBLIC m PUBLIC stdc++)
#    endif()
#endif()

# ensure preprocessor flags and module dependencies containing preprocessor flags are respected

# define the output modules directory

set_target_properties( ${PMLIB} PROPERTIES Fortran_MODULE_DIRECTORY "${ParaMonte_MOD_DIR}" )
# target_compile_definitions( ${PMLIB} INTERFACE "${FPP_FLAGS}" ) # cmake currently does not add preprocessor automatically

#string(REPLACE ";" " " FCL_FLAGS ${FCL_FLAGS})
#string(REPLACE ";" " " FC_LIB_FLAGS ${FC_LIB_FLAGS})
message( STATUS "${pmattn} ParaMonte library build compile flags: ${FCL_FLAGS} ${FC_LIB_FLAGS}" )
#target_compile_options( ${PMLIB} BEFORE PRIVATE "${FCL_FLAGS}" "${FC_LIB_FLAGS}" )
target_compile_options( ${PMLIB} BEFORE PRIVATE "${FCL_FLAGS}" "${FC_LIB_FLAGS}" )

#message( STATUS "${pmattn} ParaMonte library build linker flags: ${${FCL_FLAGS}${FL_LIB_FLAGS}}${MATLAB_LINKER_FLAGS}" )
#target_link_options( ${PMLIB} BEFORE PRIVATE ${${FCL_FLAGS}${FL_LIB_FLAGS}}${MATLAB_LINKER_FLAGS} )
message( STATUS "${pmattn} ParaMonte library build linker flags: ${${FCL_FLAGS}${FL_LIB_FLAGS}}" )
target_link_options( ${PMLIB} BEFORE PRIVATE ${${FCL_FLAGS}${FL_LIB_FLAGS}} )

install(TARGETS ${PMLIB} 
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        DESTINATION "${ParaMonte_LIB_DIR}/"
        )
