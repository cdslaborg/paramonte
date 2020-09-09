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

# test source files
set(ParaMonteTest_SOURCE_FILES
    ${ParaMonteTest_SRC_DIR}/Test_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Decoration_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Err_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_DateTime_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_String_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_System_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Path_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_FileList_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Timer_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_TimerCPU_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_RandomSeed_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_FileContents_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_File_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_CrossCorr_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Matrix_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_CorrCoef_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_TranGaus_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Math_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Misc_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_BandSpectrum_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Batse_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Statistics_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_Optimization_mod.f90
    ${ParaMonteTest_SRC_DIR}/Test_ParaMonte_mod.f90
    ${ParaMonteTest_SRC_DIR}/main.f90
)

# Define the executable in terms of the source files

# Have the .mod files placed in the lib folder
#unset(CMAKE_Fortran_MODULE_DIRECTORY)
#set(CMAKE_Fortran_MODULE_DIRECTORY "${ParaMonteTest_MOD_DIR}")

#find_library(PMLIB_PATH ${PMLIB} PATHS ${ParaMonte_OBJ_DIR})

add_executable(${ParaMonteTest_EXE_NAME} ${ParaMonteTest_SOURCE_FILES})

if (${LTYPE} MATCHES "[Ss][Tt][Aa][Tt][Ii][Cc]")
    target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC ${PMLIB})
elseif (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]")
    target_link_libraries(${ParaMonteTest_EXE_NAME} PUBLIC ${PMLIB})
endif()

set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${ParaMonteTest_MOD_DIR}" )
target_include_directories(${ParaMonteTest_EXE_NAME} PRIVATE "${ParaMonte_MOD_DIR}" PRIVATE "${ParaMonteTest_MOD_DIR}" )
set_target_properties(${ParaMonteTest_EXE_NAME} PROPERTIES LINKER_LANGUAGE "Fortran")
target_compile_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FCL_FLAGS}" )
target_link_options(${ParaMonteTest_EXE_NAME} BEFORE PRIVATE "${FCL_FLAGS}" "${FL_FLAGS}" )

# set up BLAS and LAPACK if needed

if (DEFINED BLAS AND DEFINED LAPACK)
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
#: install the test executable
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set(CMAKE_INSTALL_PREFIX "${ParaMonteTest_BIN_DIR}/")
install(TARGETS ${ParaMonteTest_EXE_NAME} 
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}"
        )

#if (${LTYPE} MATCHES "[Dd][Yy][Nn][Aa][Mm][Ii][Cc]")
#    set(CMAKE_INSTALL_PREFIX "${ParaMonteTest_BIN_DIR}/")
#    install(TARGETS ${PMLIB} 
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

set( TEST_INPUT_FILES
   "${ParaMonteTest_SRC_DIR}/input/autoCorrDataWeightedCompactResultRef.txt"
   "${ParaMonteTest_SRC_DIR}/input/autoCorrDataWeightedCompact.txt"
   "${ParaMonteTest_SRC_DIR}/input/paramonte.nml"
)
foreach( file_i ${TEST_INPUT_FILES})
    add_custom_command( TARGET ${ParaMonteTest_EXE_NAME} POST_BUILD
                        # ${CMAKE_COMMAND} -E achieves a cross-platform copy
                        COMMAND ${CMAKE_COMMAND} -E copy 
                        ${file_i}
                        # $<TARGET_FILE:${ParaMonteTest_EXE_NAME}> resolves to the full path of the library created by ParaMonteTest_EXE_NAME: "${ParaMonteTest_BIN_DIR}/", regardless of platform or compiler
                        "${ParaMonteTest_INPUT_DIR}/"
                        COMMENT "copying test data ${file_i} to the test runtime directory ..."
                        )
endforeach( file_i )
