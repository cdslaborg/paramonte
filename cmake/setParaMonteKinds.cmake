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

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# First infer the available kinds supported by the compiler. This must be done before all subsequent steps.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (DEFINED ski OR DEFINED iki OR DEFINED lki OR DEFINED cki OR DEFINED rki)
    unset(FPP_PM_KINDS)
endif()

if (${codecov_enabled} AND ${csid_is_gnu} AND NOT DEFINED rki)
    set(rki "1;2")
endif()

# Unless explicitly enforced, ensure the requested complex kinds have equivalent real kinds.
if (DEFINED cki AND NOT DEFINED rki)
    set(rki ${cki})
endif()

# Unless explicitly enforced, ensure the requested complex kinds have equivalent real kinds.
if (DEFINED rki AND NOT DEFINED cki)
    set(cki ${rki})
endif()

# This file must be executed at every cmake run because the pm_sampling@generics
# Fortran files depend on the `*_K_ENABLED` variables defined here.
# UPDATE: not anymore as cmake insanities files are now removed.
if ("${FPP_PM_KINDS}" STREQUAL "" OR NOT DEFINED FPP_PM_KINDS)

    set(MAX_COUNT_KINDS 5)

    unset(SK_COUNT)
    unset(IK_COUNT)
    unset(LK_COUNT)
    unset(CK_COUNT)
    unset(RK_COUNT)

    unset(KINDS_COUNT_LIST)

    set(KIND_SRC_DIR_PATH "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/kinds")
    set(KIND_SRC_FILE_PATH "${KIND_SRC_DIR_PATH}/getFortranCompilerKinds.F90")
    set(KIND_EXE_FILE_PATH "${KIND_SRC_DIR_PATH}/main.exe")

    message(NOTICE "${pmattn} Checking for the supported kinds by compiler: ${CMAKE_Fortran_COMPILER}")
    message(NOTICE "${pmattn} Generating file: ${KIND_SRC_FILE_PATH}")

    file(WRITE "${KIND_SRC_FILE_PATH}"
        "
        use iso_fortran_env
        write(output_unit, '(*(g0,:,'';''))', advance = 'no')   size(character_kinds), &
                                                                size(integer_kinds), &
                                                                size(logical_kinds), &
                                                                size(real_kinds)
        end
        "
        )

    try_run(PROGRAM_RUN_NOTICE
            PROGRAM_COMPILE_SUCCEEDED
            "${KIND_SRC_DIR_PATH}" # Binary dir
            "${KIND_SRC_FILE_PATH}" # source file
            COMPILE_OUTPUT_VARIABLE BUILD_OUTPUT
            RUN_OUTPUT_VARIABLE KINDS_COUNT_LIST
            )

    if (PROGRAM_RUN_NOTICE STREQUAL FAILED_TO_RUN OR NOT PROGRAM_COMPILE_SUCCEEDED)
        PrintTestCompilerStatus("Fortran" "  -- iso_fortran_env info extraction failed.")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log "The Fortran kinds extraction failed with the following output:\n${BUILD_OUTPUT}\n\n")
        string(REPLACE "\n" " \n  " CLEANED_BUILD_OUTPUT "${BUILD_OUTPUT}")
        message(FATAL_ERROR
                "${pmfatal} Failed to capture the supported kind type parameters of the compiler:\n"
                "${pmfatal} \n"
                "${pmfatal}     ${CMAKE_Fortran_COMPILER}\n"
                "${pmfatal} \n"
                "${pmfatal} Possible reasons include:\n"
                "${pmfatal} \n"
                "${pmfatal}     -#  The lack of compiler support for the `iso_fortran_env` intrinsic module.\n"
                "${pmfatal}     -#  Other invasive applications (like Dropbox, ...) have locked access to the relevant resources (particularly on Windows).\n"
                "${pmfatal}     -#  Other unknown failures in building or running executable.\n"
                "${pmfatal} \n"
                "${pmfatal} Here is the output from the compilation step:\n"
                "${pmfatal} \n"
                "${CLEANED_BUILD_OUTPUT}\n"
                "${pmfatal} \n"
                "${pmfatal} Here is the output from the running step:\n"
                "${pmfatal} \n"
                "${KINDS_COUNT_LIST}\n"
                "\n"
                )
    else()
        message(NOTICE "${pmattn} The supported kinds counts (string, integer, logical, real): ${KINDS_COUNT_LIST}")
        list(GET KINDS_COUNT_LIST 0 SK_COUNT)
        list(GET KINDS_COUNT_LIST 1 IK_COUNT)
        list(GET KINDS_COUNT_LIST 2 LK_COUNT)
        list(GET KINDS_COUNT_LIST 3 RK_COUNT)
        set(CK_COUNT ${RK_COUNT})
        message(NOTICE "${pmattn} The number of supported character kinds: ${SK_COUNT}")
        message(NOTICE "${pmattn} The number of supported integer kinds: ${IK_COUNT}")
        message(NOTICE "${pmattn} The number of supported logical kinds: ${LK_COUNT}")
        message(NOTICE "${pmattn} The number of supported complex kinds: ${CK_COUNT}")
        message(NOTICE "${pmattn} The number of supported real kinds: ${RK_COUNT}")
    endif()

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Set the user-specified kinds.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    #   The following macro defines variables of pattern `<type>K<kind>_ENABLED=1`
    #   where `kind` is an integer representing the index of desired kind value in the
    #   corresponding kinds vector of the `iso_fortran_env` intrinsic module:
    #
    #       character_kinds(:), integer_kinds(:), logical_kinds(:), real_kinds(:)
    #
    #   The placeholder <type> refers to the corresponding type:
    #
    #   `S` :   character
    #   `I` :   integer
    #   `L` :   logical
    #   `C` :   complex
    #   `R` :   real
    #
    #   All other kinds not specified remain undefined until later below in this module where
    #   they are set to zero, to be later passed as FPP macros to the compiler to build the library.

    macro(defineKindMacros)

        if (${ARGC} LESS 3)
            message(FATAL_ERROR
                    "\n"
                    "${pmfatal} Internal error occurred.\n"
                    "${pmfatal} Incorrect use of CMake macro detected.\n"
                    "${pmfatal} The number of arguments passed to the macro must be two.\n"
                    "\n"
                    )
        endif()

        set(KIND_NAME ${ARGV1})
        set(TYPE_NAME ${ARGV2})
        string(TOUPPER ${KIND_NAME} KIND_NAME)

        if ("${${ARGV0}}" MATCHES "[dD][eE][fF][aA][uU][lL][tT]" OR "${${ARGV0}}" STREQUAL "" OR NOT DEFINED ${ARGV0}) # ${ARGV0} is the name of the variable passed.

            # Define all available kinds up to maximum possible count.

            if (${${KIND_NAME}_COUNT} LESS 1)
                message(FATAL_ERROR "${pmfatal} The number of supported kinds for ${KIND_NAME} kind name is less than 1: ${${KIND_NAME}_COUNT}\n\n")
            elseif (${${KIND_NAME}_COUNT} GREATER ${MAX_COUNT_KINDS})
                message(WARNING
                        "\n"
                        "${pmwarn} The number of supported real kinds is greater than ${MAX_COUNT_KINDS}: ${${KIND_NAME}_COUNT}\n"
                        "${pmwarn} Note that the ParaMonte library is currently configured to allow \n"
                        "${pmwarn} only up to five different kinds simultaneously in the build.\n"
                        "\n"
                        )
            endif()

            unset(${ARGV0})
            foreach(COUNTER RANGE 1 ${${KIND_NAME}_COUNT})
                set(${KIND_NAME}${COUNTER}_ENABLED ${COUNTER})# CACHE STRING "The kind type parameter" FORCE)#
                if (${csid_is_gnu} AND "${KIND_NAME}" STREQUAL "SK" AND NOT "${COUNTER}" STREQUAL "1")
                    message(NOTICE "${pmwarn} Skipping kind type parameter: ${KIND_NAME}${COUNTER}_ENABLED=${${KIND_NAME}${COUNTER}_ENABLED}")
                #elseif (${csid_is_gnu} AND "${KIND_NAME}" STREQUAL "LK" AND "${COUNTER}" STREQUAL "5")
                #    message(NOTICE "${pmwarn} Skipping kind type parameter: ${KIND_NAME}${COUNTER}_ENABLED=${${KIND_NAME}${COUNTER}_ENABLED}")
                else()
                    set(FPP_PM_KINDS "${FPP_PM_KINDS}" "${KIND_NAME}${COUNTER}_ENABLED=${COUNTER}")
                    message(NOTICE "${pmattn} Enabling kind type parameter: ${KIND_NAME}${COUNTER}_ENABLED=${${KIND_NAME}${COUNTER}_ENABLED}")
                    list(APPEND ${ARGV0} ${COUNTER})
                endif()
            endforeach()

            #math(EXPR START "${${KIND_NAME}_COUNT}+1")

        elseif (NOT "${${ARGV0}}" MATCHES "[nN][oO][nN][eE]")

            # Define only the user specified kinds.

            set(KIND_INDICES_LIST "${${ARGV0}}")
            string(REPLACE " " ";" KIND_INDICES_LIST "${KIND_INDICES_LIST}")
            string(REPLACE "," ";" KIND_INDICES_LIST "${KIND_INDICES_LIST}")
            string(REPLACE "/" ";" KIND_INDICES_LIST "${KIND_INDICES_LIST}")
            while (KIND_INDICES_LIST MATCHES "  " OR KIND_INDICES_LIST MATCHES ";;" OR KIND_INDICES_LIST MATCHES ",," OR KIND_INDICES_LIST MATCHES "//")
                string(REPLACE ";;" ";" KIND_INDICES_LIST "${KIND_INDICES_LIST}") # This is needed, as REMOVE_DUPLICATES below is not enough.
            endwhile()
            list(REMOVE_DUPLICATES KIND_INDICES_LIST) # remove duplicate kind indices. This is crucial to avoid ambiguous generic interfaces at build time.

            list(LENGTH KIND_INDICES_LIST LEN_KIND_INDICES_LIST)
            if (${LEN_KIND_INDICES_LIST} LESS 1 OR ${LEN_KIND_INDICES_LIST} GREATER ${MAX_COUNT_KINDS})
                printUsage()
                message(FATAL_ERROR
                        printUsage()
                        "\n"
                        "${pmfatal} The number of specified kind type parameter indices \n"
                        "${pmfatal} cannot be less than 1 or more than ${MAX_COUNT_KINDS}.\n"
                        "\n"
                        "           ${ARGV0}=\"${${ARGV0}}\"\n"
                        "           KIND_INDICES_LIST=${KIND_INDICES_LIST}\n"
                        "           LEN_KIND_INDICES_LIST=${LEN_KIND_INDICES_LIST}\n"
                        "\n"
                        "${pmfatal} See the above guidelines for help on kind specifications.\n"
                        "\n"
                        )
            endif()

            set(AT_LEAST_ONE_KIND_DEFINED FALSE)

            foreach(COUNTER RANGE 1 ${LEN_KIND_INDICES_LIST})

                math(EXPR index "${COUNTER}-1")
                list(GET KIND_INDICES_LIST ${index} KIND_INDEX)

                if ("${KIND_INDEX}" MATCHES "^[0-9]+$" AND KIND_INDEX GREATER 0 AND KIND_INDEX LESS_EQUAL ${KIND_NAME}_COUNT)
                    set(${KIND_NAME}${COUNTER}_ENABLED ${KIND_INDEX})# CACHE STRING "The kind type parameter" FORCE)#
                    set(FPP_PM_KINDS "${FPP_PM_KINDS}" "${KIND_NAME}${COUNTER}_ENABLED=${KIND_INDEX}")
                    message(NOTICE "${pmattn} Pre-specified kind type parameters detected: ${KIND_NAME}${COUNTER}_ENABLED=${${KIND_NAME}${COUNTER}_ENABLED}")
                    set(AT_LEAST_ONE_KIND_DEFINED TRUE)
                elseif (NOT "${KIND_INDEX}" STREQUAL "")
                    message(FATAL_ERROR
                            printUsage()
                            "\n"
                            "${pmfatal} Unsupported kind type parameter index for kind ${KIND_NAME}: ${KIND_INDEX}\n"
                            "${pmfatal} The kind type parameter index must be an integer\n"
                            "${pmfatal} between 1 and the length of the constant vector of the\n"
                            "${pmfatal} corresponding type in the `iso_fortran_env` intrinsic module.\n"
                            "${pmfatal} For example, if the type of interest is `integer`, an index\n"
                            "${pmfatal} can be an integer between `1` and `size(integer_kinds)` where\n"
                            "${pmfatal} `integer_kinds` is taken from the `iso_fortran_env` intrinsic module.\n"
                            "\n"
                            )
                endif()

            endforeach()

            if (NOT AT_LEAST_ONE_KIND_DEFINED)
                message(FATAL_ERROR
                        printUsage()
                        "\n"
                        "${pmfatal} The specified kind type parameter index list ${ARGV0} must not be empty.\n"
                        "\n"
                        "   ${ARGV0}=${${ARGV0}}\n"
                        "   KIND_INDICES_LIST=${KIND_INDICES_LIST}\n"
                        "\n"
                        "${pmfatal} Read the above guidelines for further help and build directions.\n"
                        "\n"
                        )
            endif()

            #math(EXPR START "${LEN_KIND_INDICES_LIST}+1")

        endif()

        #foreach(COUNTER RANGE START ${MAX_COUNT_KINDS})
        #    set(${KIND_NAME}${COUNTER}_ENABLED ${COUNTER})
        #    message(NOTICE "Disabling kind type parameter: ${KIND_NAME}${COUNTER}_ENABLED=${${KIND_NAME}${COUNTER}_ENABLED}")
        #endforeach()

        set(${ARGV0} "${${ARGV0}}" CACHE STRING "The semicolon-separated indices of the ${TYPE_NAME} kind type parameter constant vector to include in the build." FORCE)

    endmacro()

    set(SK1_ENABLED 0)
    set(SK2_ENABLED 0)
    set(SK3_ENABLED 0)
    set(SK4_ENABLED 0)
    set(SK5_ENABLED 0)

    set(IK1_ENABLED 0)
    set(IK2_ENABLED 0)
    set(IK3_ENABLED 0)
    set(IK4_ENABLED 0)
    set(IK5_ENABLED 0)

    set(LK1_ENABLED 0)
    set(LK2_ENABLED 0)
    set(LK3_ENABLED 0)
    set(LK4_ENABLED 0)
    set(LK5_ENABLED 0)

    set(CK1_ENABLED 0)
    set(CK2_ENABLED 0)
    set(CK3_ENABLED 0)
    set(CK4_ENABLED 0)
    set(CK5_ENABLED 0)

    set(RK1_ENABLED 0)
    set(RK2_ENABLED 0)
    set(RK3_ENABLED 0)
    set(RK4_ENABLED 0)
    set(RK5_ENABLED 0)

    if (${csid_is_gnu} AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 20.0.0))
        message(NOTICE
                "${pmwarn} The GNU Fortran compiler `gfortran` is known to have bugs\n"
                "${pmwarn} with processing non-default Fortran character kinds.\n"
                "${pmwarn} As such, all non-default character kinds will be disabled.\n"
                "${pmwarn} You can override this by explicitly specifying the desired \n"
                "${pmwarn} kinds as a CMake argument like: cmake .. -Dski=\"1;2;3\""
                #"\n"
                #"${pmwarn} Additionally:\n"
                #"${pmwarn} The GNU Fortran compiler `gfortran` is known to have bugs\n"
                #"${pmwarn} with processing non-default Fortran logical kind #5.\n"
                #"${pmwarn} As such, the non-default logical kind #5 will be disabled.\n"
                #"${pmwarn} You can override this by explicitly specifying the desired \n"
                #"${pmwarn} kinds as a CMake argument like: cmake .. -Dlki=\"1;2;3;4;5\""
                )
    endif()

    unset(FPP_PM_KINDS)
    defineKindMacros(ski; SK; character)# WARNING: The argument separator MUST always be a semicolon.
    defineKindMacros(iki; IK; integer)  # WARNING: The argument separator MUST always be a semicolon.
    defineKindMacros(lki; LK; logical)  # WARNING: The argument separator MUST always be a semicolon.
    defineKindMacros(cki; CK; complex)  # WARNING: The argument separator MUST always be a semicolon.
    defineKindMacros(rki; RK; real)     # WARNING: The argument separator MUST always be a semicolon.

endif()
message(NOTICE "${pmattn} Enabled compile definitions: ${FPP_PM_KINDS}")

#   Add all enabled compile kinds to the compiler list of macro definitions.

#set(FPP_PM_KINDS "${FPP_PM_KINDS}")
#set(FPP_PM_KINDS "${FPP_PM_KINDS}" CACHE STRING "The semicolon-separated indices of the kind type parameter macros to include in the build." FORCE)
add_compile_definitions("${FPP_PM_KINDS}")