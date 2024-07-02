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
# Copy the ParaMonte library collection build scripts to the collection package subdirectories.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ("${collection}" STREQUAL "example")
    set(isCollectionExam TRUE)
else()
    set(isCollectionExam FALSE)
endif()

if ("${collection}" STREQUAL "benchmark")
    set(isCollectionBench TRUE)
else()
    set(isCollectionBench FALSE)
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The rule for the ParaMonte library ${collection} build. Build ${collection} only if requested.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# This phony target is only relevant for building and running ${collection}.
add_custom_target(${collection})

message(NOTICE "${pmattn} Setting the rules for copying the ParaMonte library ${collection} in ${lang} language to the installation directory...")
message(NOTICE "${pmattn} ${collection}BR=${${collection}BR}")
message(NOTICE "${pmattn} ${collection}PP=${${collection}PP}")

if (EXISTS "${origin}")

    if (NOT "${${collection}BR}" STREQUAL "" AND NOT "${${collection}PP}" MATCHES "[Nn][Oo][Nn][Ee]")
        set(${collection}BREnabled ON)
    else()
        set(${collection}BREnabled OFF)
    endif()

    if (${collection}BREnabled)
        message(NOTICE "${pmattn} Enabling the ParaMonte library ${collection} run as requested...")
    else()
        message(NOTICE "${pmwarn} Skipping the ParaMonte library ${collection} run as requested...")
    endif()

    # Adding a custom target runs the rule only at compile time,
    # which subsequently does not allow the addition of the build scripts to the ${collection} folders,
    # because the ${collection} folders have not been added yet. The above approach has the limitation that
    # it always requires build generation via CMake to update the ${collection} files and the build scripts.
    #add_custom_target(${collection} ALL COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Copy the ${collection} files from the ${origin} to the ${destin} directory.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    message(NOTICE  "${pmattn} Copying the paramonte::${lang} library ${collection} files to the installation directory...\n"
                    "${pmattn}      - from: ${origin}\n"
                    "${pmattn}      -   to: ${destin}")
    #add_custom_target(${collection} COMMAND ${CMAKE_COMMAND} -E copy_directory "${origin}" "${destin}")
    file(COPY "${origin}" DESTINATION "${destin}")

    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Build and run the ParaMonte ${collection}.
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if(${collection}BREnabled)

        message(NOTICE "${pmattn} Enabling the selected paramonte::${lang} library ${collection} build & run...")

        ####
        #### Check if an input list of ${collection} folders is specified.
        #### If a file path is specified, store the path in ``paramonte_${collection}BR_file``.
        ####

        unset(paramonte_${collection}BR_file)
        # Assume an absolute path specified.
        set(pathtemp "${${collection}BR}")
        if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
            set(paramonte_${collection}BR_file "${${collection}BR}")
        else()
            # Assume a relative path with respect to root directory.
            set(pathtemp "${paramonte_dir}/${${collection}BR}")
            if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                set(paramonte_${collection}BR_file "${pathtemp}")
            else()
                # Assume a relative path with respect to ${collection}/${lang} root directory.
                set(pathtemp "${paramonte_dir}/${collection}/${lang}/${${collection}BR}")
                if (EXISTS "${pathtemp}" AND NOT IS_DIRECTORY "${pathtemp}")
                    set(paramonte_${collection}BR_file "${pathtemp}")
                endif()
            endif()
        endif()

        ####
        #### Create the collection items list.
        ####

        unset(paramonte_${collection}BR_dir_list)
        if (DEFINED paramonte_${collection}BR_file)

            message(NOTICE "${pmattn} The specified ${collection} file list detected at: ${paramonte_${collection}BR_file}")
            file(STRINGS "${paramonte_${collection}BR_file}" paramonte_${collection}BR_dir_list)
            # Prepend the build directory to all supplied ${collection} folders to run.
            #list(TRANSFORM paramonte_${collection}BR_dir_list PREPEND "${destin}/")

        else()

            string(REPLACE "." ";" ${collection}List "${${collection}BR}")
            string(REPLACE "," ";" ${collection}List "${${collection}List}")
            string(REPLACE " " ";" ${collection}List "${${collection}List}")
            message(NOTICE "${pmattn} Enabling the requested paramonte::${lang} library ${collection} as outlined by the flag value \"${${collection}BR}\"")
            unset(paramonte_${collection}BR_path_list)
            file(GLOB_RECURSE paramonte_${collection}BR_path_list LIST_DIRECTORIES false RELATIVE "${origin}" "${origin}/*main${lang_ext}")

            ####
            #### Refine the list to include only those requested.
            ####

            foreach(paramonte_collection_inc_path ${paramonte_${collection}BR_path_list})
                # We don't want to keep track of the ``main.m`` in the ${collection}s root directory.
                if (NOT "${paramonte_collection_inc_path}" STREQUAL "main${lang_ext}")
                    # Remove the file name from the path and keep only the directory.
                    string(REPLACE "main${lang_ext}" "" paramonte_collection_inc_dir "${paramonte_collection_inc_path}")
                    set(paramonte_collection_inc_dir_abs "${destin}/${paramonte_collection_inc_dir}")
                    if (IS_DIRECTORY "${paramonte_collection_inc_dir_abs}")
                        foreach(bitem ${${collection}List}) # no quotations for ${collection}List!
                            string(TOLOWER "${bitem}" bitem_lower)
                            if ("${bitem_lower}" STREQUAL "all" OR "${paramonte_collection_inc_dir}" MATCHES ".*${bitem}.*")
                                list(APPEND paramonte_${collection}BR_dir_list "${paramonte_collection_inc_dir}")
                            endif()
                        endforeach()
                    else()
                        message(NOTICE "${pmwarn} The globbed ${collection} directory does not exists:")
                        message(NOTICE "${pmwarn} ")
                        message(NOTICE "${pmwarn}     paramonte_collection_inc_dir=\"${paramonte_collection_inc_dir}\"")
                        message(NOTICE "${pmwarn}     paramonte_collection_inc_path=\"${paramonte_collection_inc_path}\"")
                        message(NOTICE "${pmwarn} ")
                        message(NOTICE "${pmwarn} Something is not right. Skipping the ${collection} directory...")
                    endif()
                endif()
            endforeach()

        endif()

        ####
        #### Write the collection list to an output file for
        #### subsequent parsing by the main collection script.
        #### Write each directory path on a separate line.
        ####

        string(REPLACE ";" "\n" paramonte_${collection}BR_dir_list "${paramonte_${collection}BR_dir_list}")
        file(WRITE "${destin}/exam.list" "${paramonte_${collection}BR_dir_list}")

        ####
        #### Add the custom command to run the main collection script in the root collection directory.
        ####

        message(NOTICE "${pmattn} Enabling the selected paramonte::${lang} library ${collection} builds...")


        set(buildComment "${pmattn} Running the ${collection} post-processing script: ${Cyan} ${destin}/main${lang_ext} ${ColorReset}")
        if ("${lang}" STREQUAL "matlab")
            add_custom_command(TARGET "${collection}" POST_BUILD WORKING_DIRECTORY "${destin}" COMMAND matlab -batch main COMMENT "${buildComment}" USES_TERMINAL)
        elseif ("${lang}" STREQUAL "python")
            add_custom_command(TARGET "${collection}" POST_BUILD WORKING_DIRECTORY "${destin}" COMMAND "${Python3_EXECUTABLE}" main${lang_ext} COMMENT "${buildComment}" USES_TERMINAL)
        else()
            message(NOTICE "${pmfatal} Unrecognized unsupported language ${lang} for building the ParaMonte ${collection}.")
        endif()

    endif()

else()

    message(NOTICE
            "\n"
            "${pmwarn} The ParaMonte::${lang} ${collection} directory in the root directory of the ParaMonte library ${collection} is missing.\n"
            "${pmwarn} paramonte_collection_${lang}_dir=${origin}\n"
            "${pmwarn} Skipping the ParaMonte::${lang} library ${collection} build...\n"
            "\n"
            )

endif()