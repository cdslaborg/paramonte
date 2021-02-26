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

# test source files
set(ParaMonteTest_SOURCE_FILES
    ${ParaMonteTest_SRC_DIR}/Test_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_BandSpectrum_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Batse_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Constants_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_CorrCoef_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Cosmology_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_CrossCorr_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_DateTime_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Decoration_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Err_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_File_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_FileContents_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_FileList_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Integration_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Kmeans_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_KmeansOOP_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Knn_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Math_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Matrix_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_MultiSkewNorm_mod.f90
    #${ParaMonteTest_SRC_DIR}/Test_Partition_mod.f90
    #${ParaMonteTest_SRC_DIR}/Test_PartitionKmeans_mod.f90
    #${ParaMonteTest_SRC_DIR}/Test_PartitionVarDen_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_PartitionMaxDen_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_PartitionMinVol_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Misc_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Optimization_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDRAM_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDISE_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecBase_smod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecMCMC_smod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecDRAM_smod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_runSampler_smod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Parallelism_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDRAM_ChainFileContents_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDISE_ChainFileContents_mod.f90
    #${ParaMonteTest_SRC_DIR}/Test_ParaNest_ChainFileContents_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDRAM_RefinedChain_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaDISE_RefinedChain_mod.f90
    #${ParaMonteTest_SRC_DIR}/Test_ParaNest_RefinedChain_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Path_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_RandomSeed_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Sort_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_StarFormation_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Statistics_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_GeoCyclicFit_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_String_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_System_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Timer_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_TimerCPU_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_TranGaus_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Unique_mod.f90
    ${ParaMonteTest_SRC_DIR}/main.f90
    #${ParaMonteTest_SRC_DIR}/Test_EconomicsToolbox_mod.f90
)

set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_PartitionMaxDen_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_Partition_mod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_PartitionMinVol_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_Partition_mod.inc.f90)

set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_runSampler_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_runSampler_smod.inc.f90)

set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecBase_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecBase_smod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecMCMC_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecMCMC_smod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecDRAM_smod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_mod@Test_SpecDRAM_smod.inc.f90)

set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaDRAM_ChainFileContents_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_ChainFileContents_mod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaDISE_ChainFileContents_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_ChainFileContents_mod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaNest_ChainFileContents_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_ChainFileContents_mod.inc.f90)

set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaDRAM_RefinedChain_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_RefinedChain_mod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaDISE_RefinedChain_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_RefinedChain_mod.inc.f90)
set_source_files_properties(${ParaMonteTest_SRC_DIR}/Test_ParaNest_RefinedChain_mod.f90 PROPERTIES OBJECT_DEPENDS ${ParaMonteTest_SRC_DIR}/Test_ParaXXXX_RefinedChain_mod.inc.f90)

# Define the executable in terms of the source files

# Have the .mod files placed in the lib folder
#unset(CMAKE_Fortran_MODULE_DIRECTORY)
#set(CMAKE_Fortran_MODULE_DIRECTORY "${ParaMonteTest_MOD_DIR}")

#find_library(PMLIB_PATH ${PMLIB_NAME} PATHS ${ParaMonte_OBJ_DIR})

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#: add test target, only if testing is enabled
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (BASIC_TEST_ENABLED OR SAMPLER_TEST_ENABLED OR CODECOV_ENABLED)

    add_executable(${ParaMonteTest_EXE_NAME} ${ParaMonteTest_SOURCE_FILES})

    set(Test_SRC_BENCH ${ParaMonteTest_SRC_DIR}/___test_benchmn___.f90)
    if (EXISTS ${Test_SRC_BENCH})
        target_sources(${ParaMonteTest_EXE_NAME} PRIVATE ${Test_SRC_BENCH})
    endif()

    if (CODECOV_ENABLED AND gnu_compiler)
        target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC gcov)
    endif()

    if (${LTYPE} MATCHES "[Ss][Tt][Aa][Tt][Ii][Cc]")
        target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC ${PMLIB_NAME})
    elseif (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]" OR ${LTYPE} MATCHES "[Ss][Hh][Aa][Rr][Ee][Dd]")
        target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC ${PMLIB_NAME})
    endif()

    set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${ParaMonteTest_MOD_DIR}" )
    target_include_directories(${ParaMonteTest_EXE_NAME} PRIVATE "${ParaMonte_MOD_DIR}" PRIVATE "${ParaMonteTest_MOD_DIR}" )
    set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES LINKER_LANGUAGE "Fortran")
    target_compile_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FPP_FLAGS}" "${FC_FLAGS}" "${FCL_FLAGS}")
    target_link_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FCL_FLAGS}" "${FL_FLAGS}" )

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
                              COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                              LINK_FLAGS    "${OpenMP_Fortran_FLAGS}"
                              )
    endif()

    if (MPI_ENABLED)
        set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES
                              COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                              LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
        target_include_directories(${ParaMonteTest_EXE_NAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
        target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC ${MPI_Fortran_LIBRARIES})
    endif ()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #: report the test executable build options
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    get_target_property(PMLIB_TEST_INCLUDE_DIRECTORIES   ${ParaMonteTest_EXE_NAME} INCLUDE_DIRECTORIES)
    get_target_property(PMLIB_TEST_COMPILE_DEFINITIONS   ${ParaMonteTest_EXE_NAME} COMPILE_DEFINITIONS)
    get_target_property(PMLIB_TEST_COMPILE_OPTIONS       ${ParaMonteTest_EXE_NAME} COMPILE_OPTIONS)
    get_target_property(PMLIB_TEST_COMPILE_FLAGS         ${ParaMonteTest_EXE_NAME} COMPILE_FLAGS)

    get_target_property(PMLIB_TEST_LINK_DIRECTORIES  ${ParaMonteTest_EXE_NAME} LINK_DIRECTORIES)
    get_target_property(PMLIB_TEST_INSTALL_NAME_DIR  ${ParaMonteTest_EXE_NAME} INSTALL_NAME_DIR)
    get_target_property(PMLIB_TEST_LINK_OPTIONS      ${ParaMonteTest_EXE_NAME} LINK_OPTIONS)
    get_target_property(PMLIB_TEST_LINK_FLAGS        ${ParaMonteTest_EXE_NAME} LINK_FLAGS)

    message( STATUS "${pmattn} ParaMonte library test target include directories: ${PMLIB_TEST_INCLUDE_DIRECTORIES}" )
    message( STATUS "${pmattn} ParaMonte library test target compile definitions: ${PMLIB_TEST_COMPILE_DEFINITIONS}" )
    message( STATUS "${pmattn} ParaMonte library test target compile options    : ${PMLIB_TEST_COMPILE_OPTIONS}" )
    message( STATUS "${pmattn} ParaMonte library test target compile flags      : ${PMLIB_TEST_COMPILE_FLAGS}" )
    message( STATUS "${pmattn} ParaMonte library test target link directories   : ${PMLIB_TEST_LINK_DIRECTORIES}" )
    message( STATUS "${pmattn} ParaMonte library test target install name dir   : ${PMLIB_TEST_INSTALL_NAME_DIR}" )
    message( STATUS "${pmattn} ParaMonte library test target link options       : ${PMLIB_TEST_LINK_OPTIONS}" )
    message( STATUS "${pmattn} ParaMonte library test target link flags         : ${PMLIB_TEST_LINK_FLAGS}" )

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #: install the test executable
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    set(CMAKE_INSTALL_PREFIX "${ParaMonteTest_BIN_DIR}/")
    install(TARGETS ${ParaMonteTest_EXE_NAME} 
            CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
            RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}"
            )

    #if (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]")
    #    set(CMAKE_INSTALL_PREFIX "${ParaMonteTest_BIN_DIR}/")
    #    install(TARGETS ${PMLIB_NAME} 
    #            CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
    #            DESTINATION "${CMAKE_INSTALL_PREFIX}"
    #            )
    #endif()

    message( STATUS "${pmattn} CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )

    #if (NOT CMAKE_FORTRAN_LINK_EXECUTABLE)
    #    set(CMAKE_FORTRAN_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_FLAGS> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
    #    message( STATUS "CMAKE_FORTRAN_LINK_EXECUTABLE: ${CMAKE_FORTRAN_LINK_EXECUTABLE}")
    #endif ()

    # copy input data for testing to the bin folder

    enable_testing()

    file(GLOB TEST_INPUT_FILES "${ParaMonteTest_SRC_INPUT_DIR}/*")

    message( STATUS "${pmattn} copying the tests' input files" )
    message( STATUS "${pmattn} from: ${ParaMonteTest_SRC_INPUT_DIR}/" )
    message( STATUS "${pmattn}   to: ${ParaMonteTest_BIN_INPUT_DIR}/" )

    foreach( file_i ${TEST_INPUT_FILES})
        add_custom_command( TARGET ${ParaMonteTest_EXE_NAME} POST_BUILD
                            # ${CMAKE_COMMAND} -E achieves a cross-platform copy
                            COMMAND ${CMAKE_COMMAND} -E copy 
                            ${file_i}
                            # $<TARGET_FILE:${ParaMonteTest_EXE_NAME}> resolves to the full path of the library created by ParaMonteTest_EXE_NAME: "${ParaMonteTest_BIN_DIR}/", regardless of platform or compiler
                            "${ParaMonteTest_BIN_INPUT_DIR}/"
                            COMMENT "copying test data ${file_i} to the test runtime directory ..."
                            )
    endforeach( file_i )

endif()