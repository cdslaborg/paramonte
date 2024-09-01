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

# ParaMonte library

unset(paramonte_src_main_fortran_files)
unset(paramonte_src_main_fortran_modlist)

# sterile modules

set(sterile
    pm_array
    pm_bit
    pm_blas
    pm_control
    pm_lapack
    pm_mathConst
    pm_memory
    pm_paramonte
    pm_physUnit
    pm_search
    pm_strANSI
    pm_test
    pm_timer
)
set(paramonte_src_main_fortran_modlist "${paramonte_src_main_fortran_modlist}" "${sterile}")
foreach(modname ${sterile})
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}.F90")
endforeach()

# insanity modules

set(insanity
    pm_sampling_scio
    pm_sampling_base
    pm_sampling_mcmc
    pm_sampling_dram
    pm_sampling_dise
    pm_sampling_nest
    pm_sampling_kernel
    pm_sampling_proposal
)
set(paramonte_src_main_fortran_modlist "${paramonte_src_main_fortran_modlist}" "${insanity}")
foreach(modname ${insanity})
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}.F90")
    set_source_files_properties("${paramonte_src_fortran_main_dir}/${modname}.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/${modname}.imp.F90")
endforeach()
# STUPID CMAKE OVERWRITES OBJECT DEPENDENCIES WITH EVERY CALL TO SET_SOURCE_FILES_PROPERTIES.
get_source_file_property(old_deps "${paramonte_src_fortran_main_dir}/pm_sampling_scio.F90" OBJECT_DEPENDS)
set(all_deps "${old_deps}" "${paramonte_src_fortran_main_dir}/pm_sampling_scio.inc.F90")
set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_scio.F90" PROPERTIES OBJECT_DEPENDS "${all_deps}")
get_source_file_property(old_deps "${paramonte_src_fortran_main_dir}/pm_sampling_kernel.F90" OBJECT_DEPENDS)
set(all_deps "${old_deps}" "${paramonte_src_fortran_main_dir}/pm_sampling_kernel.inc.F90" "${paramonte_src_fortran_main_dir}/pm_sampling_kernel.inc.inc.F90")
set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_kernel.F90" PROPERTIES OBJECT_DEPENDS "${all_deps}")
get_source_file_property(srcfileDepList "${paramonte_src_fortran_main_dir}/pm_sampling_kernel.F90" OBJECT_DEPENDS)
#message(FATAL_ERROR "${srcfileDepList}")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_scio.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_scio.F90")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_base.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_scio.F90")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_mcmc.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_base.F90")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_nest.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_base.F90")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_dram.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_mcmc.F90")
#set_source_files_properties("${paramonte_src_fortran_main_dir}/pm_sampling_dise.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_sampling_dram.F90")

## insanity modules
#
#foreach(COUNTER RANGE 1 5)
#    if (RK${COUNTER}_ENABLED)
#        set(temp "${paramonte_src_fortran_main_dir}/pm_sampling@generics.rk")
#        set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${temp}${COUNTER}.F90")
#        set_source_files_properties("${temp}" PROPERTIES OBJECT_DEPENDS "${temp}.imp.F90")
#        set_source_files_properties("${temp}" PROPERTIES OBJECT_DEPENDS "${temp}.inc.F90")
#    endif()
#endforeach()

# fertile modules

set(grandpa
    pm_arrayCenter
    pm_arrayChange
    pm_arrayChoice
    pm_arrayCompact
    pm_arrayCompareLex
    pm_arrayComplement
    pm_arrayCopy
    pm_arrayFill
    pm_arrayFind
    pm_arrayInit
    pm_arrayInsert
    pm_arrayMembership
    pm_arrayMerge
    pm_arrayMinMax
    pm_arrayPad
    pm_arrayRange
    pm_arrayRank
    pm_arrayRefine
    pm_arrayRemap
    pm_arrayRemove
    pm_arrayReplace
    pm_arrayResize
    pm_arrayReverse
    pm_arraySearch
    pm_arraySelect
    pm_arrayShuffle
    pm_arraySort
    pm_arraySpace
    pm_arraySplit
    pm_arrayStrip
    pm_arrayUnique
    pm_arrayVerbose
    pm_batse
    pm_bench
    pm_clustering
    pm_complexAbs
    pm_complexCompareAll
    pm_complexCompareAny
    pm_complexCompareLex
    pm_complexDiv
    pm_complexMinMax
    pm_container
    pm_cosmicRate
    pm_cosmology
    pm_dateTime
    pm_kind
    pm_io
    pm_distanceBhat
    pm_distanceEuclid
    pm_distanceHellinger
    pm_distanceKolm
    pm_distanceMahal
    pm_distBand
    pm_distBern
    pm_distBeta
    pm_distBinom
    pm_distCosRaised
    pm_distCov
    pm_distExp
    pm_distExpGamma
    pm_distGenExpGamma
    pm_distGamma
    pm_distGenGamma
    pm_distGeom
    pm_distGeomCyclic
    pm_distKolm
    pm_distLogNorm
    pm_distLogUnif
    pm_distEggBox
    pm_distMultiNorm
    pm_distNormShell
    pm_distUnifEll
    pm_distUnifPar
    pm_distNegExp
    pm_distNorm
    pm_distPareto
    pm_distPiwiPoweto
    pm_distPower
    pm_distPoweto
    pm_distPois
    pm_distUnif
    pm_distUnifSphere
    pm_io
    pm_err
    pm_except
    pm_fftnr
    pm_fftpack
    pm_ellipsoid
    pm_polation
    pm_knn
    pm_logicalCompare
    pm_math1mexp
    pm_mathBeta
    pm_mathCumPropExp
    pm_mathCumSum
    pm_mathCompare
    pm_mathDivMul
    pm_mathErf
    pm_mathExp
    pm_mathFactorial
    pm_mathFactoring
    pm_mathFisher
    pm_mathGamma
    pm_mathGammaAM
    pm_mathGammaGil
    pm_mathGammaNR
    pm_mathLog1p
    pm_mathLogAddExp
    pm_mathLogSubExp
    pm_mathLogSumExp
    pm_mathMinMax
    pm_mathNumSys
    pm_mathRoot
    pm_mathRootTest
    pm_mathRound
    pm_mathSqrt
    pm_mathSubAdd
    pm_mathSum
    pm_matrixMulTri
    pm_mathUnsigned
    pm_matrixChol
    pm_matrixClass
    pm_matrixCopy
    pm_matrixDet
    pm_matrixIndex
    pm_matrixInit
    pm_matrixInv
    pm_matrixLUP
    pm_matrixMulAdd
    pm_matrixPack
    pm_matrixSubset
    pm_matrixTrace
    pm_matrixTrans
    pm_matrixUpdate
    pm_optimization
    pm_option
    pm_os
    pm_parallelism
    pm_polynomial
   #pm_procPois
    pm_quadPack
    pm_quadRomb
    pm_quadTest
    pm_sampleAffinity
    pm_sampleACT
    pm_sampleCCF
    pm_sampleCor
    pm_sampleCov
    pm_sampleECDF
    pm_sampleMean
    pm_sampleNorm
    pm_sampleQuan
    pm_sampleScale
    pm_sampleShift
    pm_sampleVar
    pm_sampleWeight
    pm_sampling
    pm_statest
    pm_str
    pm_strASCII
    pm_swap
    pm_sysInfo
    pm_sysPath
    pm_sysShell
    pm_val2complex
    pm_val2int
    pm_val2logical
    pm_val2real
    pm_val2str
    pm_ziggurat
    )
set(paramonte_src_main_fortran_modlist "${paramonte_src_main_fortran_modlist}" "${grandpa}")
foreach(modname ${grandpa})
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}.F90")
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}@routines.F90")
    set(incsrc "${paramonte_src_fortran_main_dir}/${modname}@routines.inc.F90")
    if (EXISTS "${incsrc}")
        set_source_files_properties("${paramonte_src_fortran_main_dir}/${modname}@routines.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/${modname}@routines.inc.F90")
    endif()
endforeach()

# bizarre cases.

set(resizes
    pm_arrayRebill
    pm_arrayRefill
    pm_arrayRebind
    )
set(paramonte_src_main_fortran_modlist "${paramonte_src_main_fortran_modlist}" "${resizes}")
foreach(modname ${resizes})
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}.F90")
    set(paramonte_src_main_fortran_files "${paramonte_src_main_fortran_files}" "${modname}@routines.F90")
    set_source_files_properties("${paramonte_src_fortran_main_dir}/${modname}@routines.F90" PROPERTIES OBJECT_DEPENDS "${paramonte_src_fortran_main_dir}/pm_arrayResize@routines.inc.F90")
endforeach()

#    ${paramonte_src_fortran_main_dir}/EconomicsToolbox_pmod.F90
#    ${paramonte_src_fortran_main_dir}/pm_procPois.F90
#    ${paramonte_src_fortran_main_dir}/pm_histogram.F90
#    ${paramonte_src_fortran_main_dir}/pm_distMultiNormGrid.F90
#    ${paramonte_src_fortran_main_dir}/pm_distNormShellMix.F90
#    ${paramonte_src_fortran_main_dir}/DistMultiNormGenExpGammaMix_pmod.F90
#    ${paramonte_src_fortran_main_dir}/pm_distMultiSkewNorm.F90
#    ${paramonte_src_fortran_main_dir}/ClusteredPoint_pmod.F90
#    ${paramonte_src_fortran_main_dir}/pm_partition.F90
#    ${paramonte_src_fortran_main_dir}/pm_partition@MaxDen_smod.F90
#    ${paramonte_src_fortran_main_dir}/pm_partition@MinVol_smod.F90
#    #${paramonte_src_fortran_main_dir}/PartitionKmeans_pmod.F90
#    #${paramonte_src_fortran_main_dir}/PartitionVarDen_pmod.F90
#    #${paramonte_src_fortran_main_dir}/pm_partitionRecursiveMaxDen.F90
#    #${paramonte_src_fortran_main_dir}/pm_partitionRecursiveMinVol.F90
#    #${paramonte_src_fortran_main_dir}/pm_partitionMaxDen.F90
#    ${paramonte_src_fortran_main_dir}/pm_tranGaus.F90

# Clean up the source list to add to the target.

set(depscan_enabled ON)
if ("${mod}" STREQUAL "reset")
    set(mod "${mod_old}" CACHE STRING "default module list to compile." FORCE)
elseif (DEFINED paramonte_src_main_fortran_paths_refined)
    if ("${mod}" STREQUAL "" OR "${mod}" STREQUAL "${mod_old}") # if not set ever or if equal to the past value.
        set(depscan_enabled OFF)
    endif()
endif()

if (depscan_enabled AND "${mod}" STREQUAL "")
    if (NOT "${lang}" STREQUAL "fortran")
        set(mod "pm_sampling")
    endif()
endif()

set(mod_old "${mod}" CACHE INTERNAL "copy of user-specified module lists to compile" FORCE)

# Generate the requested module source list.

if (depscan_enabled)

    if ("${mod}" STREQUAL "" OR "${mod}" STREQUAL "all")

        set(paramonte_src_main_fortran_files_refined "${paramonte_src_main_fortran_files}")

    else()

        message(NOTICE "${pmattn} ParaMonte Fortran library source module file search pattern: ${mod}")

        unset(paramonte_src_main_fortran_files_refined)
        foreach(pattern ${mod})
            foreach(srcfile ${paramonte_src_main_fortran_files})
                if ("${srcfile}" MATCHES ".*${pattern}.*")
                    set(paramonte_src_main_fortran_files_refined "${paramonte_src_main_fortran_files_refined}" "${srcfile}")
                endif()
            endforeach()
        endforeach()

        #### Search all files for dependencies inclusion.

        set(depCheckFileList "${paramonte_src_main_fortran_files_refined}")
        list(LENGTH depCheckFileList depCheckFileListLen)

        while (NOT "${depCheckFileListLen}" EQUAL 0)

            message(NOTICE "${pmattn} depCheckFileListLen, depCheckFileList: ${depCheckFileListLen}, ${depCheckFileList}")

            #### depCheckFileListNew contains the list of all newly added sources file whose contents must be checked.
            #### As such, depCheckFileListNew is always a subset of paramonte_src_main_fortran_files_refined at any moment at runtime.

            unset(depCheckFileListNew)
            foreach(srcfile ${depCheckFileList})

                setLineList("${paramonte_src_fortran_main_dir}/${srcfile}" lineList)

                #### Append the dependency source lines.

                get_source_file_property(srcfileDepList "${paramonte_src_fortran_main_dir}/${srcfile}" OBJECT_DEPENDS)
                foreach(srcdep ${srcfileDepList})
                    if (EXISTS "${srcdep}")
                        message(NOTICE "${pmattn} additional object dependency file to scan for ${srcfile}: ${srcfileDepList}")
                        setLineList("${srcdep}" depLineList)
                        set(lineList "${lineList}" "${depLineList}")
                    endif()
                endforeach()

                #### Now a list of all source lines for the current source file ${srcfile} is in the list `lineList`.
                #### Search for module name mentions in each line of the file.

                #### Comparing lower-case may be faster than case-insensitive pattern matching. conclusion: It is not.
                #string(REPLACE " " "" lineList "${lineList}")
                #string(TOLOWER "${lineList}" lineList)

                foreach(line ${lineList})

                    #if ("${line}" MATCHES ".*usepm_.*")
                    if ("${line}" MATCHES ".*[Uu][Ss][Ee][ ]*[Pp][Mm]_.*")

                        #### Strip blanks at the beginning and the end.

                        string(STRIP "${line}" line)

                        #### Exclude comment lines.
                        #### \todo
                        #### The following approach comes with the caveat that we are also
                        #### including modules names that may appear in strings or end of line comments.
                        #### This should be fixed in future.

                        #### Find the comment symbol in the current line.
                        #### If line is not a full comment (i.e., the zeroth element is not `!`), we do a full search for module `use` presence.

                        string(FIND "${line}" "!" icol)
                        if(NOT "${icol}" EQUAL 0)

                            #### Not a comment line.

                            string(TOLOWER "${line}" line)
                            foreach(modname ${paramonte_src_main_fortran_modlist})

                                string(TOLOWER "${modname}" modnameLower)
                                if ("${line}" MATCHES ".*${modnameLower}.*")

                                    #### include the dependency source file if it is not already included.

                                    if (NOT "${modname}.F90" IN_LIST paramonte_src_main_fortran_files_refined)

                                        set(paramonte_src_main_fortran_files_refined "${paramonte_src_main_fortran_files_refined}" "${modname}.F90")
                                        set(depCheckFileListNew "${depCheckFileListNew}" "${modname}.F90")

                                        #### Include any submodules of the module as well.
                                        #### Here, we assume submodule source file naming always follows the pattern: modname@*.F90

                                        foreach(item ${paramonte_src_main_fortran_files})
                                            if ("${item}" MATCHES "${modname}@.*")
                                                set(paramonte_src_main_fortran_files_refined "${paramonte_src_main_fortran_files_refined}" "${item}")
                                                set(depCheckFileListNew "${depCheckFileListNew}" "${item}")
                                            endif()
                                        endforeach()

                                    endif()

                                endif()

                            endforeach()

                        endif()

                    endif()

                endforeach()

            endforeach()

            set(depCheckFileList "${depCheckFileListNew}")
            list(LENGTH depCheckFileList depCheckFileListLen)

        endwhile()

    endif()

    #### add source directory to file names.

    set(paramonte_src_main_fortran_paths_refined "${paramonte_src_main_fortran_files_refined}" CACHE INTERNAL "source file list" FORCE)

endif()

message(NOTICE "${pmattn} ParaMonte Fortran library source files to compile: ${paramonte_src_main_fortran_paths_refined}")
set(paramonte_src_main_fortran_paths_refined_abs ${paramonte_src_main_fortran_paths_refined})
list(TRANSFORM paramonte_src_main_fortran_paths_refined_abs PREPEND "${paramonte_src_fortran_main_dir}/")
#message(FATAL_ERROR "${pmattn} ParaMonte Fortran library source files to compile: ${paramonte_src_main_fortran_paths_refined_abs}")

#unset(paramonte_src_main_fortran_paths_refined_abs)
#foreach(srcfile ${paramonte_src_main_fortran_files})
#    foreach(srcfile_selected ${paramonte_src_main_fortran_paths_refined})
#        if ("${srcfile}" STREQUAL "${srcfile_selected}")
#            set(paramonte_src_main_fortran_paths_refined_abs "${paramonte_src_main_fortran_paths_refined_abs}" "${paramonte_src_fortran_main_dir}/${srcfile}")
#            break()
#        endif()
#    endforeach()
#endforeach()

# enable debug info on IO source files

#if (${csid_is_gnu})
#    set(SRC_LEVEL_DEBUG_FLAG "-g")
#elseif(${csid_is_intel})
#    if (WIN32)
#        set_source_files_properties(${paramonte_src_fortran_main_dir}/pm_sampling@Input_smod.F90 PROPERTIES COMPILE_FLAGS "/debug")
#        set_source_files_properties(${paramonte_src_fortran_main_dir}/pm_sampling@Input_smod.F90 PROPERTIES COMPILE_FLAGS "/debug")
#    elseif (UNIX OR APPLE)
#        set_source_files_properties(${paramonte_src_fortran_main_dir}/pm_sampling@Input_smod.F90 PROPERTIES COMPILE_FLAGS "-g")
#        set_source_files_properties(${paramonte_src_fortran_main_dir}/pm_sampling@Input_smod.F90 PROPERTIES COMPILE_FLAGS "-g") # ${INTEL_Fortran_DEBUG_FLAGS}
#    endif()
#endif()

# search for static libraries

function(find_static_library libname OUT)

    if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    elseif (UNIX)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()

    find_library(
        FOUND_${libname}_STATIC
        "${libname}"
        )

    if (FOUND_${libname}_STATIC)
        get_filename_component(ABS_FILE ${FOUND_${libname}_STATIC} ABSOLUTE)
        set(${OUT} ${ABS_FILE} PARENT_SCOPE)
    else()
        message(NOTICE
                "\n"
                "${pmwarn}\n"
                "${pmwarn} Unable to find ${libname}${CMAKE_FIND_LIBRARY_SUFFIXES}\n"
                "${pmwarn} Please report any failures to shahmoradi@utexas.edu "
                "\n"
                )
        set(${OUT} "" PARENT_SCOPE)
    endif()

endfunction()

unset(FL_FLAGS_PMLIB)
unset(FCL_FLAGS_PMLIB)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add OpenBLAS if it exists.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#target_include_directories(${target} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set IPO/LTO flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

unset(FCL_FLAGS_IPO)
if ("${build}" STREQUAL "ipo" OR ${build} STREQUAL "tuned" OR ${build} STREQUAL "native")

    if(${csid_is_gnu})

        set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
        -fopt-info-all=GFortranOptReport.txt
        -finline-functions      # consider all functions for inlining, even if they are not declared inline.
       #-fwhole-program         # allow the compiler to make assumptions on the visibility of the symbols leading to more aggressive optimization decisions.
        -flto=3                 # enable interprocedural optimization between files in parallel on 3 processors.
        )

    elseif(${csid_is_intel})

        if (WIN32)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            /Qopt-report:2          # generate optimization report. Level 2 is the default. Use 5 for the greatest details and 0 for no report.
            /Qipo                   # enable interprocedural optimization between files.
            /Qip                    # determines whether additional interprocedural optimizations for single-file compilation are enabled.
            )

        elseif (UNIX)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            -qopt-report=2          # generate optimization report. Level 2 is the default. Use 5 for the greatest details and 0 for no report.
            -inline-forceinline     # Instructs the compiler to force inlining of functions suggested for inlining whenever the compiler is capable doing so.
            -finline-functions      # enables function inlining for single file compilation.
            -ipo                    # enable interprocedural optimization between files.
            -ip                     # determines whether additional interprocedural optimizations for single-file compilation are enabled.
            )

        endif()

    else()

        # Optional IPO. Do not use IPO if it's not supported by compiler.
        check_ipo_supported(RESULT result OUTPUT output)
        if(result)
          set_property(TARGET "${libname}" PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
        else()
          message(FATAL_ERROR "IPO is not supported: ${output}")
        endif()

    endif()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set the native vs. tuned multiple dispatch optimization flags
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ("${build}" STREQUAL "tuned")

    if(${csid_is_gnu})

        set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
        -mtune=intel # Produce code optimized for the most current Intel processors (Haswell and Silvermont as of GCC version 12).
        )

    elseif(${csid_is_intel})

        if (WIN32)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            /QaxCOMMON-AVX512,CORE-AVX512,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2,HASWELL,KNL,SKYLAKE,SKYLAKE-AVX512,BROADWELL,SILVERMONT
            )

        elseif (APPLE)

            #   SSE3,SSE2
            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            -axCOMMON-AVX512,CORE-AVX512,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,HASWELL,SKYLAKE,SKYLAKE-AVX512,BROADWELL
            )

        elseif (UNIX)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            -axCOMMON-AVX512,CORE-AVX512,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2,HASWELL,KNL,SKYLAKE,SKYLAKE-AVX512,BROADWELL
            )

        endif()

    else()

        # Warn the user about the lack of tuned support for generic compilers.
        message(WARNING
                "${pmwarn} Unrecognized unsupported compiler detected."
                "Please contact the ParaMonte developers for assistance with tuned builds via unsupported compilers.")

    endif()

elseif ("${build}" STREQUAL "native")

    if(${csid_is_gnu})

        set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
        -march=native # Produce code optimized for the most current Intel processors (Haswell and Silvermont as of GCC version 12).
        )

    elseif(${csid_is_intel})

        if (WIN32)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            /QxHost
            )

        elseif (APPLE)

            #   SSE3,SSE2
            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            -xHost
            )

        elseif (UNIX)

            set(FCL_FLAGS_IPO "${FCL_FLAGS_IPO}"
            -xHost
            )

        endif()

    else()

        # Warn the user about the lack of tuned support for generic compilers.
        message(WARNING
                "${pmwarn} Unrecognized unsupported compiler detected."
                "Please contact the ParaMonte developers for assistance with tuned builds via unsupported compilers.")

    endif()

endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add the C library first.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# For C source files, we will create a separate Object Library (https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/Object-Library)

set(libnamec "${libname}c")
unset(paramonte_src_main_c_paths)
set(paramonte_src_main_c_paths "${paramonte_src_c_main_dir}/pm_sys.c")
add_library(${libnamec} OBJECT "${paramonte_src_main_c_paths}")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add the library and compile/link flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


if (DLL_ENABLED)

    # The following does not seem to work.
    #if (${csid_is_gnu} AND UNIX AND NOT APPLE)
    #    #set(CMAKE_SKIP_RPATH ON CACHE BOOL "Skip adding a build-dir rpath")
    #    set(CMAKE_BUILD_WITH_INSTALL_RPATH ON CACHE BOOL "Build using the install rpath")
    #    set(CMAKE_BUILD_RPATH_USE_ORIGIN ON CACHE BOOL "Use relative rpaths")
    #    set(CMAKE_INSTALL_RPATH "$ORIGIN")
    #endif()

    add_library("${libname}" SHARED "${paramonte_src_main_fortran_paths_refined_abs}" $<TARGET_OBJECTS:${libnamec}>)

    #if (${csid_is_gnu})
        #set_target_properties( "${libname}" PROPERTIES COMPILE_FLAGS "-static-libgfortran")
        #set_target_properties( "${libname}" PROPERTIES LINK_FLAGS "-static-libgfortran")
        #set_target_properties( "${libname}" PROPERTIES COMPILE_FLAGS "-Wl,--whole-archive -lgfortran -static -lgfortran")
        #set(CMAKE_SHARED_LINKER_FLAGS "-static-libgfortran")
        #set_target_properties( "${libname}" PROPERTIES LINK_FLAGS "-static-libgfortran")
        #if (APPLE)
            #set_target_properties( "${libname}" PROPERTIES LINK_FLAGS "-Wl,--whole-archive")
        #elseif(UNIX)
            #set_target_properties( "${libname}" PROPERTIES LINK_FLAGS "-Wl,-static-libgfortran")
        #endif()

        #set_target_properties( "${libname}" PROPERTIES SKIP_BUILD_RPATH FALSE) # don't skip the full RPATH for the build tree
        #set_target_properties( "${libname}" PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE) # when building, don't use the install RPATH for the build
        #set_target_properties( "${libname}" PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE) # if TRUE, append the automatically determined parts of RPATH
        #set_target_properties( "${libname}" PROPERTIES INSTALL_RPATH "./") # the RPATH to be used when installing the library
        #set_target_properties( "${libname}" PROPERTIES LINK_FLAGS "-Wl,-rpath,'$ORIGIN'") # if uncommented, this will override all of the above!
    #endif()

   #if (${csid_is_gnu})
   #    find_static_library(libgfortran LIBGFORTRAN_PATH)
   #    find_static_library(libquadmath LIBQUADMATH_PATH)
   #    if (${build} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries("${libname}" debug ${LIBGFORTRAN_PATH})
   #        target_link_libraries("${libname}" debug ${LIBQUADMATH_PATH})
   #    else() #if (${build} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries("${libname}" optimized ${LIBGFORTRAN_PATH})
   #        target_link_libraries("${libname}" optimized ${LIBQUADMATH_PATH})
   #    endif()
   #elseif(${csid_is_intel})
   #    find_static_library(libifport LIBIFPORT_PATH)
   #    if (${build} MATCHES "[dD][eE][bB][uU][gG]")
   #        target_link_libraries("${libname}" debug ${LIBIFPORT_PATH})
   #    else() #if (${build} MATCHES "[rR][eE][lL][eE][aA][sS][eE]")
   #        target_link_libraries("${libname}" optimized ${LIBIFPORT_PATH})
   #    endif()
   #endif()

    set_target_properties("${libname}" PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
    set_target_properties("${libname}" PROPERTIES BUILD_RPATH_USE_ORIGIN ON)

    unset(rpath_prop)
    if (APPLE)
        set(rpath_prop "@loader_path")
    elseif(UNIX)
        set(rpath_prop "$ORIGIN")
    endif()

    if (MPI_ENABLED)
        set(rpath_prop "${rpath_prop};${MPI_Fortran_INCLUDE_DIRS}")
    endif()

    if (DEFINED rpath_prop)
        set_target_properties("${libname}" PROPERTIES INSTALL_RPATH "${rpath_prop}")
    endif()

    #### set up the library-specific compile link flags

    if(${csid_is_gnu})

        set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
        -shared
        -fPIC
        )

        set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
        -shared
        )

        #if (DEPS_ENABLED)
            #set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            #-static-libgfortran
            #-static-libgcc
            #-static
            #)
        #endif()

    elseif(${csid_is_intel})

        if (WIN32)

            set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
            #/threads " # these flags are actually included by default in recent ifort implementations
            /libs:dll
            )

            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            /static # On Windows* systems, option /static is equivalent to option /MT.
            /dll
            )

        elseif(UNIX)

            set(FCL_FLAGS_PMLIB "${FCL_FLAGS_PMLIB}"
            -fpic # Request compiler to generate position-independent code.
            )

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

            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            -qopenmp-link=static
            -static-intel
            -static
            )

            #if (DEPS_ENABLED)
                set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
                -static-libstdc++
                -static-libgcc
                )
            #endif()

        endif()

        #### Include MKL

        #if (CMAKE_BUILD_TYPE MATCHES "Release|RELEASE|release")
        #    #find_package(MKL REQUIRED)
        #    target_link_directories("${libname}" PUBLIC ${MKL_LIBRARIES})
        #    target_include_directories("${libname}" PUBLIC ${MKL_INCLUDE_DIRS})
        #    target_link_libraries("${libname}" PUBLIC mkl_intel_lp64 mkl_sequential mkl_core)
        #endif()

    else()

        option(BUILD_SHARED_LIBS ON)

    endif()

else() # static lib

    add_library("${libname}" STATIC "${paramonte_src_main_fortran_paths_refined_abs}" $<TARGET_OBJECTS:${libnamec}>)
    if (${csid_is_intel})
        if (WIN32)
            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            -static
            )
        else()
            set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
            -qopenmp-link=static
            -static-intel
            -static
            )
            #if (DEPS_ENABLED)
                set(FL_FLAGS_PMLIB "${FL_FLAGS_PMLIB}"
                -static-libstdc++
                -static-libgcc
                )
            #endif()
        endif()
    endif()

endif()

#### Bypass the (potential) GNU gfortran 13.2 bug for PPC arch in `pm_kind`:
#foreach(srcfile ${paramonte_src_main_fortran_paths_refined_abs})
#    if ("${srcfile}" MATCHES ".*pm_kind.*")
#        set_source_files_properties("${srcfile}" PROPERTIES COMPILE_FLAGS "-fno-range-check")
#    endif()
#endforeach()

#set_target_properties("${libname}" PROPERTIES VERSION "${ParaMonteVersion}" SOVERSION "${ParaMonteVersionMajor}")

if(OpenMP_Fortran_FOUND)
    target_link_libraries("${libname}" PUBLIC OpenMP::OpenMP_Fortran)
endif()

if (OpenBLAS_ENABLED)
    # Warning: The following openblas library name will have to updated to `openblas_64`
    # if the compiler options for long integer kinds for array size are enabled.
    target_link_libraries("${libname}" PUBLIC openblas)
    #set_target_properties("${libname}" PROPERTIES LINK_FLAGS "-Wl,-rpath,.")
elseif(BLAS_ENABLED)
endif()

if (MPI_ENABLED)
    #include_directories(${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries("${libname}" PUBLIC ${MPI_Fortran_LIBRARIES})
    target_include_directories("${libname}" PUBLIC ${MPI_Fortran_INCLUDE_PATH})
endif()

if (MATLAB_ENABLED AND Matlab_FOUND)
    message(NOTICE "${pmattn} Matlab_INCLUDE_DIRS=${Matlab_INCLUDE_DIRS}")
    target_include_directories("${libname}" PUBLIC "${Matlab_INCLUDE_DIRS}")
    target_link_libraries("${libname}" PUBLIC "${Matlab_MEX_LIBRARY}" "${Matlab_MX_LIBRARY}")# "${Matlab_LIBRARIES}"
    message(NOTICE "${pmattn} Matlab_MEX_LIBRARY=${Matlab_MEX_LIBRARY}")
    message(NOTICE "${pmattn} Matlab_MX_LIBRARY=${Matlab_MX_LIBRARY}")
    #if (NOT CAF_ENABLED AND NOT MPI_ENABLED)
        #allow matlab console print only in serial MATLAB applications.
        target_compile_definitions("${libname}" PUBLIC "MEXPRINT_ENABLED=1")
        #add_compile_definitions("MEXPRINT_ENABLED=1")
    #endif()
else()
    target_compile_definitions("${libname}" PUBLIC "MEXPRINT_ENABLED=0")
    #add_compile_definitions("MEXPRINT_ENABLED=0")
endif()

#if (NOT WIN32)
#    if (DEFINED MATLAB_ROOT_DIR)
#        target_include_directories("${libname}" PUBLIC "${MATLAB_ROOT_DIR}/extern/include")
#        target_link_directories("${libname}" BEFORE
#        PUBLIC "$(MATLAB_ROOT_DIR)/bin/${MATLAB_PLATFORM_SUBDIR}"
#        PUBLIC "$(MATLAB_ROOT_DIR)/extern/bin/${MATLAB_PLATFORM_SUBDIR}"
#        )
#        target_link_libraries("${libname}" PUBLIC mx PUBLIC mex PUBLIC mat PUBLIC m PUBLIC stdc++)
#    endif()
#endif()

# ensure preprocessor flags and module dependencies containing preprocessor flags are respected

# define the output modules directory

set_target_properties("${libname}" PROPERTIES Fortran_MODULE_DIRECTORY "${paramonte_bld_inc_dir}")
# target_compile_definitions( "${libname}" INTERFACE "${FPP_FLAGS}") # cmake currently does not add preprocessor automatically

#string(REPLACE ";" " " FCL_FLAGS ${FCL_FLAGS})
#message(NOTICE "${pmattn} ParaMonte library build flags: ${FPP_FLAGS} ${FCL_FLAGS}")
target_compile_options("${libname}" BEFORE PRIVATE "${FC_FLAGS}" "${FCL_FLAGS}" "${FCL_FLAGS_PMLIB}" "${FCL_FLAGS_IPO}" "${fcf}")

#message(NOTICE "${pmattn} ParaMonte library build linker flags: ${${FCL_FLAGS}${FL_FLAGS_PMLIB}}${MATLAB_LINKER_FLAGS}")
#target_link_options( "${libname}" BEFORE PRIVATE ${${FCL_FLAGS}${FL_FLAGS_PMLIB}}${MATLAB_LINKER_FLAGS} )
#message(NOTICE "${pmattn} ParaMonte library build linker flags: ${FCL_FLAGS} ${FL_FLAGS_PMLIB}")
#target_link_options( "${libname}" BEFORE PRIVATE ${${FCL_FLAGS}${FL_FLAGS}${FL_FLAGS_PMLIB}} )
target_link_options("${libname}" BEFORE PRIVATE "${FCL_FLAGS}" "${FL_FLAGS}" "${FL_FLAGS_PMLIB}" "${FCL_FLAGS_IPO}" "${flf}")

get_target_property(PMLIB_INCLUDE_DIRECTORIES   "${libname}" INCLUDE_DIRECTORIES)
get_target_property(PMLIB_COMPILE_DEFINITIONS   "${libname}" COMPILE_DEFINITIONS)
get_target_property(PMLIB_COMPILE_OPTIONS       "${libname}" COMPILE_OPTIONS)
get_target_property(PMLIB_COMPILE_FLAGS         "${libname}" COMPILE_FLAGS)

get_target_property(PMLIB_LINK_DIRECTORIES  "${libname}" LINK_DIRECTORIES)
get_target_property(PMLIB_INSTALL_NAME_DIR  "${libname}" INSTALL_NAME_DIR)
get_target_property(PMLIB_LINK_OPTIONS      "${libname}" LINK_OPTIONS)
get_target_property(PMLIB_LINK_FLAGS        "${libname}" LINK_FLAGS)

message(NOTICE "${pmattn} target include directories: ${PMLIB_INCLUDE_DIRECTORIES}")
message(NOTICE "${pmattn} target compile definitions: ${PMLIB_COMPILE_DEFINITIONS}")
message(NOTICE "${pmattn} target compile options    : ${PMLIB_COMPILE_OPTIONS}")
message(NOTICE "${pmattn} target compile flags      : ${PMLIB_COMPILE_FLAGS}")
message(NOTICE "${pmattn} target link directories   : ${PMLIB_LINK_DIRECTORIES}")
message(NOTICE "${pmattn} target install name dir   : ${PMLIB_INSTALL_NAME_DIR}")
message(NOTICE "${pmattn} target link options       : ${PMLIB_LINK_OPTIONS}")
message(NOTICE "${pmattn} target link flags         : ${PMLIB_LINK_FLAGS}")
message(NOTICE "${pmattn} CMAKE INSTALL PATH        : ${paramonte_bld_lib_dir}")

if (${codecov_enabled} AND ${csid_is_gnu})
    target_link_libraries("${libname}" PUBLIC gcov)
elseif (${codecov_enabled})
    message(FATAL_ERROR "${pmfatal} Code Coverage is currently only supported with the GNU compilers.\n")
endif()

set(SRC_BENCHM ${paramonte_src_fortran_main_dir}/___benchmn___.F90)
if (EXISTS "${SRC_BENCHM}")
    find_package(LAPACK REQUIRED)
    set(BLAS_FOUND TRUE CACHE INTERNAL "BLAS was found" FORCE)
    set(LAPACK_FOUND TRUE CACHE INTERNAL "LAPACK was found" FORCE)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "BLAS LIBS" FORCE)
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "LAPACK LIBS" FORCE)
    target_sources("${libname}" PRIVATE "${SRC_BENCHM}")
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Ensure C source flags are set properly.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set_source_files_properties("${paramonte_src_main_c_paths}" PROPERTIES LINK_FLAGS "")
set_source_files_properties("${paramonte_src_main_c_paths}" PROPERTIES LINK_OPTIONS "")
set_source_files_properties("${paramonte_src_main_c_paths}" PROPERTIES COMPILE_FLAGS "")
set_source_files_properties("${paramonte_src_main_c_paths}" PROPERTIES COMPILE_OPTIONS "")
set_source_files_properties("${paramonte_src_main_c_paths}" PROPERTIES COMPILE_DEFINITIONS "")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Install the library in the lib folder of the build directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

install(TARGETS "${libname}"
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        RUNTIME_DEPENDENCY_SET pmdeps
        )
#install(RUNTIME_DEPENDENCY_SET pmdeps
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
#        )

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Recipe for installing the runtime dependencies.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (DLL_ENABLED AND DEPS_ENABLED)
    install(CODE [[
        ####
        #### For intel compilers, we take an aggressive approach to find the dependencies
        #### as several libraries in the Intel OneAPI collection do not link statically.
        ####
        if ("${CMAKE_INSTALL_PREFIX}" MATCHES ".*[Ii][Nn][Tt][Ee][Ll].*")
            execute_process(COMMAND which ifx OUTPUT_VARIABLE ifcpath)
            string(REGEX REPLACE "\n$" "" ifcpath "${ifcpath}")
            if (IS_DIRECTORY "${ifcpath}" OR NOT EXISTS "${ifcpath}")
                execute_process(COMMAND which ifx OUTPUT_VARIABLE ifcpath)
                string(REGEX REPLACE "\n$" "" ifcpath "${ifcpath}")
                if (IS_DIRECTORY "${ifcpath}" OR NOT EXISTS "${ifcpath}")
                    unset(ifcpath)
                endif()
            endif()
        endif()
        set(ilibdir ".")
        if (DEFINED ifcpath)
            get_filename_component(ilibdir "${ifcpath}" DIRECTORY)
            get_filename_component(ilibdir "${ilibdir}" DIRECTORY)
            set(ilibdir "${ilibdir}/lib")
            if (NOT IS_DIRECTORY "${ilibdir}")
                message(NOTICE "${pmwarn} Failed to find the Intel compilers library path from: ${ifcpath}")
                message(NOTICE "${pmwarn} The expected Intel compilers library path is: ${ilibdir}")
                set(ilibdir ".")
            endif()
        endif()
        ####
        #### Now recursively find dependencies.
        ####
        function(install_library_with_deps LIBRARY)
            file(INSTALL
                DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
                TYPE SHARED_LIBRARY
                FOLLOW_SYMLINK_CHAIN
                FILES "${LIBRARY}"
            )
            file(GET_RUNTIME_DEPENDENCIES
                LIBRARIES ${LIBRARY}
                RESOLVED_DEPENDENCIES_VAR resolved_deps
                UNRESOLVED_DEPENDENCIES_VAR unresolved_deps
                CONFLICTING_DEPENDENCIES_PREFIX colflicting
                POST_EXCLUDE_REGEXES ".*MATLAB.*" ".*matlab.*"
                DIRECTORIES "${ilibdir}"
            )
            foreach(FILE ${resolved_deps})
                if(NOT IS_SYMLINK ${FILE})
                    install_library_with_deps(${FILE})
                endif()
            endforeach()
            foreach(FILE ${unresolved_deps})
                message(NOTICE "${pmwarn} Unresolved dependency from ${LIBRARY}: ${FILE}")
            endforeach()
            foreach(FILE ${colflicting_FILENAMES})
                message(NOTICE "${pmwarn} Conflicting dependency: ${FILE}")
            endforeach()
        endfunction()
        file(GET_RUNTIME_DEPENDENCIES
            LIBRARIES $<TARGET_SONAME_FILE:libparamonte>
            CONFLICTING_DEPENDENCIES_PREFIX colflicting
            UNRESOLVED_DEPENDENCIES_VAR unresolved_deps
            RESOLVED_DEPENDENCIES_VAR resolved_deps
            POST_EXCLUDE_REGEXES ".*MATLAB.*" ".*matlab.*"
            DIRECTORIES "${ilibdir}"
        )
        foreach(FILE ${resolved_deps})
            install_library_with_deps(${FILE})
        endforeach()
        foreach(FILE ${unresolved_deps})
            message(NOTICE "${pmwarn} Unresolved dependency: ${FILE}")
        endforeach()
        foreach(FILE ${colflicting_FILENAMES})
            message(NOTICE "${pmwarn} Conflicting dependency: ${FILE}")
        endforeach()
    ]])
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Add the executable for build-time FPP source generation.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# This is used in setParaMonteFPP.cmake and must be set right here!
if (FPP_ENABLED)
    add_executable(genfpp.exe "${paramonte_dir}/auxil/genfpp.F90")
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(genfpp.exe PUBLIC OpenMP::OpenMP_Fortran)
    endif()
    if (${csid_is_gnu})
        #target_compile_options(genfpp.exe BEFORE PRIVATE "${FC_FLAGS}" "${FCL_FLAGS}")
        target_compile_options(genfpp.exe BEFORE PRIVATE "-std=gnu" "-ffree-line-length-none")
    endif()
    include(setParaMonteFPP)
endif()

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## Copy the auxiliary files to the package directory. Build examples only if requested.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
##   LICENSE
#
#set(origin "${paramonte_dir}/LICENSE.md")
#set(destin "${paramonte_bld_pkg_dir}/LICENSE.md")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library LICENSE.md file to the installation directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyFileLICENSE ALL COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
#add_dependencies("${libname}" copyFileLICENSE)
#
##   README
#
#set(origin "${paramonte_src_${lang}_dir}/README.md")
#set(destin "${paramonte_bld_pkg_dir}/README.md")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library README.md file to the installation directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyFileREADME ALL COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
#add_dependencies("${libname}" copyFileREADME)
#
##   CHANGES
#
#set(origin "${paramonte_dir}/CHANGES.md")
#set(destin "${paramonte_bld_pkg_dir}/CHANGES.md")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library CHANGES.md file to the installation directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyFileCHANGES ALL COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
#add_dependencies("${libname}" copyFileCHANGES)
#
##   VERSION
#
#set(origin "${paramonte_dir}/VERSION.md")
#set(destin "${paramonte_bld_pkg_dir}/auxil/VERSION.md")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library VERSION.md file to the installation directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyFileVERSION ALL COMMAND ${CMAKE_COMMAND} -E copy "${origin}" "${destin}")
#add_dependencies("${libname}" copyFileVERSION)
#
##   INCLUDE
#
#if ("${lang}" STREQUAL "c" OR "${lang}" STREQUAL "cpp")
#
#    set(origin "${paramonte_src_${lang}_dir}/paramonte.*")
#    set(destin "${paramonte_bld_pkg_dir}/")
#    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the installation directory...\n"
#                    "${pmattn}      - from: ${origin}\n"
#                    "${pmattn}      -   to: ${destin}")
#    add_custom_target(copyDirInclude ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#    add_dependencies("${libname}" copyDirInclude)
#
#elseif ("${lang}" STREQUAL "fortran")
#
#    set(origin "${paramonte_bld_inc_dir}")
#    set(destin "${paramonte_bld_pkg_dir}/inc")
#    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the installation directory...\n"
#                    "${pmattn}      - from: ${origin}\n"
#                    "${pmattn}      -   to: ${destin}")
#    add_custom_target(copyDirInclude ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#    add_dependencies("${libname}" copyDirInclude)
#
#elseif ("${lang}" STREQUAL "julia" OR "${lang}" STREQUAL "matlab" OR "${lang}" STREQUAL "python" OR "${lang}" STREQUAL "r")
#
#    # Copy library files
#
#    set(origin "${paramonte_src_${lang}_dir}/paramonte")
#    set(destin "${paramonte_bld_pkg_dir}/paramonte")
#    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library include files to the installation directory...\n"
#                    "${pmattn}      - from: ${origin}\n"
#                    "${pmattn}      -   to: ${destin}")
#    add_custom_target(copyDirInclude ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#    add_dependencies("${libname}" copyDirInclude)
#
#    # Copy auxiliary files
#
#    set(origin "${paramonte_src_${lang}_dir}/auxil")
#    set(destin "${paramonte_bld_pkg_dir}/auxil")
#    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library auxiliary files to the installation directory...\n"
#                    "${pmattn}      - from: ${origin}\n"
#                    "${pmattn}      -   to: ${destin}")
#    add_custom_target(copyDirAuxil ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#    add_dependencies("${libname}" copyDirAuxil)
#
#    if ("${lang}" STREQUAL "python")
#
#        set(origin "${paramonte_src_${lang}_dir}/setup")
#        set(destin "${paramonte_bld_pkg_dir}/setup")
#        message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library setup files to the installation directory...\n"
#                        "${pmattn}      - from: ${origin}\n"
#                        "${pmattn}      -   to: ${destin}")
#        add_custom_target(copyDirSetup ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#        add_dependencies("${libname}" copyDirSetup)
#
#    endif()
#
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copy the ParaMonte library examples to the package directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#message(NOTICE "${pmattn} Setting the rule for copying the ParaMonte library examples in ${lang} language to the installation directory...")
#set(origin "${paramonte_example_dir}/${lang}")
#set(destin "${paramonte_bld_pkg_dir}/example")
#if (EXISTS "${origin}")
#
#    # Copy the example files.
#
#    message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library example files to the installation directory...\n"
#                    "${pmattn}      - from: ${origin}\n"
#                    "${pmattn}      -   to: ${destin}")
#    add_custom_target(copyDirExample ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#    #add_dependencies("${libname}" copyDirExample)
#
#    # Loop over all directories to add the build scripts.
#
#    if ("${build}" STREQUAL "debug")
#        set(paramonte_example_build_type "debug")
#    else()
#        set(paramonte_example_build_type "release")
#    endif()
#
#    if (${MPI_ENABLED})
#        set(paramonte_example_par_type ".mpi")
#    else()
#        unset(paramonte_example_par_type)
#    endif()
#
#    set(paramonte_example_build_script_bash  "build.${csid}.${paramonte_example_build_type}${paramonte_example_par_type}.sh")
#    set(paramonte_example_build_script_batch "build.${csid}.${paramonte_example_build_type}${paramonte_example_par_type}.bat")
#
#    # Loop over all directories to add the build scripts.
#
#    if ("${lang}" STREQUAL "fortran")
#        if (csid_is_gnu OR csid_is_intel)
#            setSubDirList(subDirList "${destin}")
#            set(counter 1)
#    #message(FATAL_ERROR "${pmattn} destin: ${destin}")
#    message(FATAL_ERROR "${pmattn} subDirList: ${subDirList}")
#    message(FATAL_ERROR "${pmattn} Bash build script: ${paramonte_example_build_script_bash}")
#    message(NOTICE "${pmattn} Batch build script: ${paramonte_example_build_script_batch}")
#            foreach(subDir ${subDirList})
#                setSubDirList(subSubDirList "${destin}/${subDir}")
#                foreach(subSubDir ${subSubDirList})
#                    set(example_current_dir "${destin}/${subDir}/${subSubDir}")
#                    message(NOTICE "${pmattn} Adding the build scripts to example at: ${example_current_dir}")
#                    #file(COPY "${paramonte_example_dir}/auxil/${lang}/${paramonte_example_build_script_bash}" DESTINATION "${example_current_dir}"
#                    #FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#                    #file(RENAME "${example_current_dir}/${paramonte_example_build_script_bash}" "${example_current_dir}/build.sh")
#                    add_custom_target(copyBuildBash${counter} ALL COMMAND ${CMAKE_COMMAND} -E copy
#                    "${paramonte_example_dir}/auxil/${lang}/${paramonte_example_build_script_bash}" "${example_current_dir}/build.sh")
#                    add_dependencies("${libname}" copyBuildBash${counter})
#                    math(EXPR counter "${counter}+1")
#                    if (WIN32)
#                        #file(COPY "${paramonte_example_dir}/auxil/${lang}/${paramonte_example_build_script_batch}" DESTINATION "${example_current_dir}"
#                        #FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#                        #file(RENAME "${example_current_dir}/${paramonte_example_build_script_batch}" "${example_current_dir}/build.bat")
#                        add_custom_target(copyBuildBatch${counter} ALL COMMAND ${CMAKE_COMMAND} -E copy
#                        "${paramonte_example_dir}/auxil/${lang}/${paramonte_example_build_script_batch}" "${example_current_dir}/build.bat")
#                        add_dependencies("${libname}" copyBuildBatch${counter})
#                        math(EXPR counter "${counter}+1")
#                    endif()
#                endforeach()
#            endforeach()
#        else()
#            message(WARNING "${pmwarn} Unrecognized unsupported compiler detected. ParaMonte examples build scripts will likely not function.")
#        endif()
#    else()
#        message(FATAL_ERROR "${pmfatal} Internal CMake error. Unsupported example language.")
#    endif()
#
#    ## Set the first line of the build scripts.
#    #
#    #set(example_build_precompile_bash "#!/bin/bash")
#    #set(example_build_precompile_batch "set PATH=..\\..\\..\\lib;%PATH%")
#    #
#    ## Set the library name for example build scripts.
#    #
#    #set(example_build_compile_bash "main.F90 ../../../lib/libparamonte_*")
#    #if (WIN32 AND "${lib}" STREQUAL "shared" AND ${csid_is_intel})
#    #    set(example_build_compile_batch "main.F90 ..\\..\\..\\lib\\libparamonte_*.lib")
#    #else()
#    #    set(example_build_compile_batch "main.F90 ..\\..\\..\\lib\\libparamonte_*")
#    #endif()
#    #
#    ## Generate the generic build scripts.
#    #
#    #if ("${lang}" STREQUAL "fortran")
#    #
#    #    if (${csid_is_gnu})
#    #
#    #        set(example_build_compile_bash "-cpp -ffree-line-length-none -Wl,-rpath,../../../lib -I../../../inc ${example_build_compile_bash}")
#    #        if ("${build}" STREQUAL "debug")
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O0 -g -fcheck=all -fbacktrace ${example_build_compile_bash}")
#    #        elseif ("${build}" STREQUAL "testing")
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O2 ${example_build_compile_bash}")
#    #        else()
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O3 ${example_build_compile_bash}")
#    #        endif()
#    #        set(example_build_compile_bash "${example_build_compile_bash} -o main.exe")
#    #        set(example_build_compile_batch "${example_build_compile_bash}")
#    #
#    #    elseif (${csid_is_intel})
#    #
#    #        set(example_build_compile_bash "-standard-semantics -fpp -Wl,-rpath,../../../lib -I../../../inc  ${example_build_compile_bash}")
#    #        if ("${build}" STREQUAL "debug")
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O0 -g3 -CB -debug full -traceback -check all -fpe0 ${example_build_compile_bash}")
#    #        elseif ("${build}" STREQUAL "testing")
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O2 ${example_build_compile_bash}")
#    #        else()
#    #            set(example_build_compile_bash "${CMAKE_Fortran_COMPILER} -O3 ${example_build_compile_bash}")
#    #        endif()
#    #        set(example_build_compile_bash "${example_build_compile_bash} -o main.exe")
#    #
#    #        set(example_build_compile_batch "/standard-semantics /fpp /I:..\\..\\..\\include  ${example_build_compile_batch}")
#    #        if ("${build}" STREQUAL "debug")
#    #            set(example_build_compile_batch "${CMAKE_Fortran_COMPILER} /Od /debug:full /CB /Qinit:snan,arrays /warn:all /gen-interfaces /traceback /check:all /fpe-all:0 /Qtrapuv ${example_build_compile_batch}")
#    #        elseif ("${build}" STREQUAL "testing")
#    #            set(example_build_compile_batch "${CMAKE_Fortran_COMPILER} /O2 ${example_build_compile_batch}")
#    #        else()
#    #            set(example_build_compile_batch "${CMAKE_Fortran_COMPILER} /O3 ${example_build_compile_batch}")
#    #        endif()
#    #        set(example_build_compile_batch "${example_build_compile_batch} /exe:main.exe")
#    #
#    #    else()
#    #
#    #        message(WARNING "${pmwarn} Unrecognized unsupported compiler detected. ParaMonte examples build scripts will likely not function.")
#    #
#    #    endif()
#    #
#    #else()
#    #
#    #    message(FATAL_ERROR "${pmfatal} Internal CMake error. Unsupported example language.")
#    #
#    #endif()
#    #
#    ## Set the last line of the build scripts.
#    #
#    #if (${MPI_ENABLED})
#    #    set(example_build_postcompile_bash "mpiexec -n 3 ./main.exe")
#    #    if (${csid_is_intel})
#    #        set(example_build_postcompile_batch "mpiexec -localonly -n 3 main.exe")
#    #    else()
#    #        set(example_build_postcompile_batch "mpiexec -n 3 main.exe")
#    #    endif()
#    #else()
#    #    set(example_build_postcompile_bash "./main.exe")
#    #    set(example_build_postcompile_batch "main.exe")
#    #endif()
#    #
#    ## Generate the full example build scripts.
#    #
#    #set(example_build_bash "${example_build_precompile_bash}\n${example_build_compile_bash}\n${example_build_postcompile_bash}")
#    #set(example_build_batch "${example_build_precompile_batch}\n${example_build_compile_batch}\n${example_build_postcompile_batch}")
#    #file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/build.bat" "${example_build_batch}")
#    #file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/build.sh" "${example_build_bash}")
#    #message(NOTICE "${pmattn} Bash build script:")
#    #message(NOTICE "${example_build_bash}")
#    #message(NOTICE "${pmattn} Batch build script:")
#    #message(NOTICE "${example_build_batch}")
#    #
#    ## Loop over all directories to add the build scripts.
#    #
#    #setSubDirList(subDirList "${destin}")
#    #foreach(subDir ${subDirList})
#    #    setSubDirList(subSubDirList "${destin}/${subDir}")
#    #    set(example_current_dir "../../../../lib/libparamonte*.so")
#    #    foreach(subSubDir ${subSubDirList})
#    #        set(example_current_dir "${destin}/${subDir}/${subSubDir}")
#    #        message(NOTICE "${pmattn} Adding the build scripts to example at: ${example_current_dir}")
#    #        #file(WRITE "${example_current_dir}/build.sh" "${example_build_bash}")
#    #        file(COPY "${CMAKE_CURRENT_BINARY_DIR}/build.sh" DESTINATION "${example_current_dir}/"
#    #        FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#    #        #file(CHMOD "${example_current_dir}/build.sh" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#    #        if (WIN32)
#    #            #file(WRITE "${example_current_dir}/build.bat" "${example_build_batch}")
#    #            file(COPY "${CMAKE_CURRENT_BINARY_DIR}/build.bat" DESTINATION "${example_current_dir}/"
#    #            FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#    #            #file(CHMOD "${example_current_dir}/build.bat" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#    #        endif()
#    #    endforeach()
#    #endforeach()
#
#endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copy the binary files from the build package directory to the binary installation directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#if ("${ddir}" STREQUAL "")
#    set(ddir "${CMAKE_SOURCE_DIR}/bin")
#    message(NOTICE "${pmattn} Unspecified ParaMonte binary install directory `ddir`.")
#    message(NOTICE "${pmattn} Defining and generating the ParaMonte build directory: ddir=\"${ddir}\"")
#else()
#    message(NOTICE "${pmattn} User-specified ParaMonte binary install directory detected. ddir=\"${ddir}\"")
#    # Ensure the path is absolute.
#    if (NOT IS_ABSOLUTE "${ddir}")
#        set(ddir "${CMAKE_SOURCE_DIR}/${ddir}")
#        message(NOTICE "${pmattn} Absolute path for the user-specified ParaMonte binary install directory: ddir=\"${ddir}\"")
#    endif()
#endif()
##set(ENV{ddir} "${ddir}")
#
#set(paramonte_bin_current_dir "${ddir}/${libname}")
#set(ENV{paramonte_bin_current_dir} "${paramonte_bin_current_dir}")
#message(NOTICE "${pmattn} All binary files will be installed at paramonte_bin_current_dir=\"${paramonte_bin_current_dir}\"")
#if (EXISTS "${paramonte_bin_current_dir}/")
#    message(NOTICE "${pmattn} Removing the existing paramonte_bin_current_dir=\"${paramonte_bin_current_dir}\"")
#    file(REMOVE_RECURSE "${paramonte_bin_current_dir}/")
#else()
#    message(NOTICE "${pmattn} Generating the non-existing paramonte_bin_current_dir=\"${paramonte_bin_current_dir}\"")
#    file(MAKE_DIRECTORY "${paramonte_bin_current_dir}")
#endif()
#
##install(TARGETS "${libname}"
##        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
##        DESTINATION "${paramonte_bin_current_dir}/lib"
##        )
#
#set(origin "${paramonte_bld_lib_dir}")
#set(destin "${paramonte_bld_pkg_dir}/lib")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library files to the package directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyDirLib ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#add_dependencies("${libname}" copyDirLib)
#
#set(origin "${paramonte_bld_pkg_dir}")
#set(destin "${paramonte_bin_current_dir}")
#message(NOTICE  "${pmattn} Setting the rule for copying the paramonte::${lang} library package files to the installation directory...\n"
#                "${pmattn}      - from: ${origin}\n"
#                "${pmattn}      -   to: ${destin}")
#add_custom_target(copyDirDeploy ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
#add_dependencies("${libname}" copyDirDeploy)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The source for the ParaMonte library examples. Build examples only if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#if (${csid_is_gnu} OR ${csid_is_intel})
#
#    if ("${exam}" STREQUAL "none" OR NOT DEFINED exam)
#
#        message(NOTICE "${pmattn} Skipping the ParaMonte library examples build as requested...")
#
#    else()
#
#        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        # Infer the shell compatibility for building ParaMonte examples
#        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#        if (DEFINED shell)
#            string(TOLOWER "${shell}" shell)
#            if ("${shell}" STREQUAL "posix")
#                set(paramonte_example_build_script "build.sh")
#                message(NOTICE "${pmattn} Assuming POSIX-compliant runtime shell as requested for build examples...")
#            elseif ("${shell}" STREQUAL "windows")
#                set(paramonte_example_build_script "build.bat")
#                message(NOTICE "${pmattn} Assuming Windows-compliant runtime shell as requested for build examples...")
#            else()
#                message(FATAL_ERROR "${pmfatal} Unrecognized runtime shell compliance: shell=${shell}\n"
#                                    "${pmfatal} It is highly unusual to fail the CMake scripts as done here.\n"
#                                    "${pmfatal} The ParaMonte CMake scripts are likely compromised.\n"
#                                    "${pmfatal} Please report this issue to the library developers.")
#            endif()
#        else()
#            set(paramonte_example_build_script "build.sh")
#            message(WARNING "${pmwarn} Unspecified runtime shell compliance.\n"
#                            "${pmwarn} This can happen only if you are building the library by invoking CMake directly.\n"
#                            "${pmwarn} The CMake scripts will assume POSIX-compliant runtime shell to build the examples.\n"
#                            "${pmwarn} Proceeding with no guarantee of examples build success...")
#        endif()
#
#        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#        # Build the ParaMonte examples
#        #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#        message(NOTICE "${pmattn} Enabling the selected paramonte::${lang} library example build...")
#        unset(paramonte_example_file)
#        if (EXISTS "${exam}")
#            set(paramonte_example_file "${exam}")
#        else()
#            set(paramonte_example_file "${paramonte_dir}/${exam}")
#            #file(TO_NATIVE_PATH "${paramonte_example_file}" paramonte_example_file)
#            if (NOT EXISTS "${paramonte_example_file}")
#                unset(paramonte_example_file)
#            endif()
#        endif()
#
#        if (DEFINED paramonte_example_file)
#            file(STRINGS "${paramonte_example_file}" paramonte_example_rel_dir_list)
#            set(counter 1)
#            foreach(paramonte_example_rel_dir ${paramonte_example_rel_dir_list})
#                set(paramonte_example_abs_dir "${paramonte_bld_pkg_dir}/example/${paramonte_example_rel_dir}")
#                set(paramonte_example_build_script_path "${paramonte_example_abs_dir}/${paramonte_example_build_script}")
#                message(NOTICE "${pmattn} Enabling the example build at: ${paramonte_example_abs_dir}")
#                if (EXISTS "${paramonte_example_build_script_path}" AND NOT IS_DIRECTORY "${paramonte_example_build_script_path}")
#                    if ("${shell}" STREQUAL "posix")
#                        add_custom_target(buildExample${counter} ALL COMMAND sh "-c" "cd ${paramonte_example_abs_dir} && ${paramonte_example_build_script_path}")
#                    else()
#                        add_custom_target(buildExample${counter} ALL COMMAND CMD "/C" "cd ${paramonte_example_abs_dir} && ${paramonte_example_build_script_path}")
#                    endif()
#                    add_dependencies("${libname}" buildExample${counter})
#                    math(EXPR counter "${counter}+1")
#                endif()
#            endforeach()
#        else()
#            printUsage()
#            message(FATAL_ERROR
#                    "\n"
#                    "${pmfatal} Unsupported user-specified ParaMonte library examples. exam=${exam}\n"
#                    "${pmfatal} Follow the guidelines above to appropriately specify the testing mode or drop the option.\n"
#                    "\n"
#                    )
#        endif()
#
#    endif()
#
#endif()

#message(NOTICE "CMAKE_Fortran_LINK_EXECUTABLE: ${CMAKE_Fortran_LINK_EXECUTABLE}")
