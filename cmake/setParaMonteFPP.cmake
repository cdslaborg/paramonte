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

#   The following cleans up FPP source files for inclusion in the final package directory.

if (${csid_is_gnu} OR ${csid_is_intel})

    # Damn gfortran changes its behavior for setting the FPP file extension between V. 10 and newer versions.
    # Therefore, the following approach does not work.
    # Instead, we will consider all possible file extensions and search for the corresponding files at runtime.
    #if (${csid_is_gnu})
    #    set(fppFileExtensionOld ".f90" ".fii" ".fi")
    #elseif (${csid_is_intel})
    #    set(fppFileExtensionOld ".i90")
    #endif()
    set(fppFileExtensionOldList ".i90" ".fii" ".fi" ".f90" ".F90" ".F90.f90" ".I90")

    # Make the FPP package directory.

    message(NOTICE "${pmattn} paramonte_bld_pkg_fpp_dir: ${paramonte_bld_pkg_fpp_dir}" )
    if (NOT EXISTS "${paramonte_bld_pkg_fpp_dir}/")
        file(MAKE_DIRECTORY "${paramonte_bld_pkg_fpp_dir}")
    elseif ("${fresh}" MATCHES ".*${fpp}.*")
        file(REMOVE_RECURSE "${paramonte_bld_pkg_fpp_dir}/")
        file(MAKE_DIRECTORY "${paramonte_bld_pkg_fpp_dir}")
    endif()

    set(fppFileExtensionNew ".F90")
    set(fppRootDir "${paramonte_bld_obj_dir}")
    message(NOTICE "${pmattn} fppRootDir=${fppRootDir}")
    message(NOTICE "${pmattn} fppFileExtensionNew=${fppFileExtensionNew}")

    ####
    #### Find the correct fpp file extension.
    ####

    unset(fppFileDirOld)
    unset(fppFileExtensionOld)
    foreach(fppext ${fppFileExtensionOldList})
        #### The unique paramonte module file name to search for.
        #### No other source file should contain pm_kind.
        set(modname "pm_kind")
        message(NOTICE "${pmattn} checking for availability of FPP file extension: \"${fppext}\"")
        message(NOTICE "${pmattn} Globbing search pattern: \"${fppRootDir}/*${modname}*${fppext}\"")
        file(GLOB_RECURSE srclist "${fppRootDir}/*${modname}*${fppext}")
        message(NOTICE "${pmattn} Globbing results: \"${srclist}\"")
        foreach(srcpath ${srclist})
            if (EXISTS "${srcpath}")
                get_filename_component(fppFileDirOld "${srcpath}" DIRECTORY)
                set(fppFileExtensionOld "${fppext}")
                break()
            endif()
        endforeach()
        if (DEFINED fppFileExtensionOld)
            break()
        endif()
    endforeach()

    if ("${fppFileDirOld}" STREQUAL "" OR "${fppFileExtensionOld}" STREQUAL "")

        message(NOTICE "${pmwarn} Failed to identify the FPP source file extension by searching in: ${fppRootDir}")
        message(NOTICE "${pmwarn} Skipping the FPP preprocessed source file generation...")

    else()

        message(NOTICE "${pmattn} fppFileDirOld=${fppFileDirOld}")
        message(NOTICE "${pmattn} fppFileExtensionOld=${fppFileExtensionOld}")
        #message(FATAL_ERROR "")

        #### This phony target is only relevant for building and running the code for cleaning the FPP source files.

        add_custom_target(cleanupFPP)
        add_dependencies(package cleanupFPP) # ensure cleanupFPP builds before package.

        unset(fpplist)
        if (TRUE)
            file(GLOB fpplist "${fppFileDirOld}/*${fppFileExtensionOld}")
        else()
            if (TRUE)
                set(fpplist "${paramonte_src_main_fortran_files_refined}")
                #### This approach requires this file to be included within
                #### the CMakelist file of paramonte_src_fortran_main_dir directory.
                #add_dependencies(cleanupFPP "${libname}")
                #add_dependencies(package cleanupFPP) # ensure cleanupFPP builds before package.
            else()
                #add_dependencies(package cleanupFPP)
                # First create a list of file names to copy.
                file(GLOB srclist "${paramonte_src_fortran_main_dir}/*.F90")
                ##message(NOTICE "${pmwarn} ${srclist}")
                #### Refine the list to only module and submodule source files.
                unset(fpplist)
                foreach(item ${srclist})
                    # each item is full source path.
                    if (NOT "${item}" MATCHES ".*\.inc\..*" AND NOT "${item}" MATCHES ".*\.imp\..*")
                        set(fpplist "${fpplist}" "${item}")
                    endif()
                endforeach()
            endif()
        endif()

        set(ifile 1)
        #unset(fppFileExtensionOld)
        foreach(fppfile ${fpplist})
            #if (NOT DEFINED fppFileExtensionOld)
            #    foreach(fppExtOld ${fppFileExtensionOldList})
            #        if (EXISTS "${fppFileDirOld}/${fppfile}${fppExtOld}")
            #            set(fppFileExtensionOld "${fppExtOld}")
            #            break()
            #        endif()
            #    endforeach()
            #    if (NOT DEFINED fppFileExtensionOld)
            #        message(FATAL_ERROR "${pmfatal} Failed to identify the FPP source file extension among the list of possibilities: ${fppFileExtensionOldList}")
            #    endif()
            #endif()
            #################################
            #string(LENGTH "${fppfile}" lenitem)
            #   \todo
            #   Stupid CMake does not offer an humane way of replacing substring at the end of string.
            #   Instead, we have to rely on the following method.
            #   But this is error prone as there may be multiple instances of ".F90" in the file path, not just at the end.
            #string(REPLACE "$.F90" "${fppFileExtensionOld}" filepath "${fppfile}")
            #if ("${fppfile}" MATCHES ".*.F90" OR "${fppfile}" MATCHES ".*.f90")
            if (EXISTS "${fppfile}")
                set(srcpath "${fppfile}")
                get_filename_component(fppname "${fppfile}" NAME_WE) # without extension.
                #set(srcpath "${fppFileDirOld}/${fppname}${fppFileExtensionOld}")
                set(dstpath "${paramonte_bld_pkg_fpp_dir}/${fppname}${fppFileExtensionNew}")
                message(NOTICE "${pmattn} Adding the FPP build scripts for:")
                message(NOTICE "${pmattn}     source file:${Cyan} ${srcpath} ${ColorReset}")
                message(NOTICE "${pmattn}     destin file:${Cyan} ${dstpath} ${ColorReset}")
                message(NOTICE "${pmattn}     ${fppFileDirOld}/genfpp.exe ${srcpath} ${dstpath}")
                #add_custom_target(cleanupFPP${ifile} COMMAND "${fppFileDirOld}/genfpp.exe ${srcpath} ${dstpath}")
                #add_dependencies(cleanupFPP cleanupFPP${ifile})
                add_custom_command( TARGET "cleanupFPP" POST_BUILD
                                    WORKING_DIRECTORY "${fppFileDirOld}"
                                    COMMAND genfpp.exe ${srcpath} ${dstpath}
                                    COMMENT "${pmattn} Running the FPP script on: ${Cyan} ${srcpath} ${ColorReset}"
                                    USES_TERMINAL
                                    )
                math(EXPR ifile "${ifile}+1")
            endif()
        endforeach()

    endif()

else()

    message(WARNING "${pmwarn} FPP preprocessing of sources files is currently unsupported for compiler suites other than GNU and Intel. Skipping...")

endif()