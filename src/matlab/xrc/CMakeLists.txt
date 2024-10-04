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

unset(XRC_FILES)
set(XRC_FILES
    pm_sampling
    pm_parallelism
)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set mex compiler flags flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

unset(CFLAGS)
unset(LDFLAGS)
unset(MEX_FLAGS)
unset(rpath_prop)
#set(MEX_FLAGS "-v -nojvm")
#set(CFLAGS "-fPIC -shared")
set(LDFLAGS "-fPIC -shared")

if (${csid_is_gnu} OR ${csid_is_intel})
    if (APPLE)
        set(rpath_prop "@loader_path")
        set(LDFLAGS "${LDFLAGS}" "-Wl,-rpath,@rpath -Wl,-rpath,@loader_path")
    elseif(UNIX)
        set(rpath_prop "$ORIGIN")
        set(LDFLAGS "${LDFLAGS}" "-Wl,-rpath,\$ORIGIN")
    endif()
endif()

if (${csid_is_gnu})
    set(LDFLAGS "${LDFLAGS}" "-static-libstdc++")
elseif (${csid_is_intel})
    if (APPLE OR UNIX)
        set(LDFLAGS "${LDFLAGS}" "-static-intel")
    endif()
endif()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Set shared library Fortran linker flags.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# This phony target is only relevant for building and running ${collection}.
add_custom_target(matlabmex)

unset(xtargets)
foreach(xtarget ${XRC_FILES})
    set(xlibname "${xtarget}")
    set(xrcpath "${CMAKE_CURRENT_SOURCE_DIR}/${xtarget}.c")
    #set(xlibname "${xtarget}_${build}_${mem}_${parname}_${checking}")
    matlab_add_mex(NAME "${xtarget}" SHARED SRC "${xrcpath}" OUTPUT_NAME "${xlibname}" LINK_TO "${libname}" R2018a)
    target_compile_definitions("${xtarget}" PUBLIC "OMP_ENABLED=${OMP_ENABLED}")
    set_target_properties("${xtarget}" PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
    set_target_properties("${xtarget}" PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set_target_properties("${xtarget}" PROPERTIES BUILD_RPATH_USE_ORIGIN ON)
    if (DEFINED rpath_prop)
        set_target_properties("${xtarget}" PROPERTIES INSTALL_RPATH "${rpath_prop}")
    endif()
    if (DEFINED MEX_FLAGS)
        #set_target_properties("${xtarget}" PROPERTIES COMPILE_FLAGS "${MEX_FLAGS}")
        target_compile_options("${xtarget}" PUBLIC "${MEX_FLAGS}")
    endif()
    if (DEFINED LDFLAGS)
        #set_target_properties("${xtarget}" PROPERTIES LINK_FLAGS "${LDFLAGS}")
        target_link_options("${xtarget}" PUBLIC "${LDFLAGS}")
        #get_target_property(_compile_flags "${xtarget}" COMPILE_FLAGS)
        #list(REMOVE_ITEM _compile_flags "-MD")
        #list(REMOVE_ITEM _compile_flags "-MT")
        #set_target_properties("${xtarget}" PROPERTIES COMPILE_FLAGS "${_compile_flags}")
    endif()
    set(xtargets "${xtargets}" "${xtarget}")
    add_dependencies(matlabmex "${xtarget}")
    get_target_property(MXLIB_INCLUDE_DIRECTORIES   "${xtarget}" INCLUDE_DIRECTORIES)
    get_target_property(MXLIB_COMPILE_DEFINITIONS   "${xtarget}" COMPILE_DEFINITIONS)
    get_target_property(MXLIB_COMPILE_OPTIONS       "${xtarget}" COMPILE_OPTIONS)
    get_target_property(MXLIB_COMPILE_FLAGS         "${xtarget}" COMPILE_FLAGS)
    get_target_property(MXLIB_LINK_DIRECTORIES      "${xtarget}" LINK_DIRECTORIES)
    get_target_property(MXLIB_INSTALL_NAME_DIR      "${xtarget}" INSTALL_NAME_DIR)
    get_target_property(MXLIB_LINK_OPTIONS          "${xtarget}" LINK_OPTIONS)
    get_target_property(MXLIB_LINK_FLAGS            "${xtarget}" LINK_FLAGS)
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target include directories: ${MXLIB_INCLUDE_DIRECTORIES}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target compile definitions: ${MXLIB_COMPILE_DEFINITIONS}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target compile options    : ${MXLIB_COMPILE_OPTIONS}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target compile flags      : ${MXLIB_COMPILE_FLAGS}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target link directories   : ${MXLIB_LINK_DIRECTORIES}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target install name dir   : ${MXLIB_INSTALL_NAME_DIR}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target link options       : ${MXLIB_LINK_OPTIONS}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex target link flags         : ${MXLIB_LINK_FLAGS}${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}MATLAB Mex CMAKE INSTALL PATH        : ${CMAKE_INSTALL_PREFIX}/lib${ColorReset}")
    message(NOTICE "${pmattn} ${BoldYellow}Matlab_VERSION                       : ${Matlab_VERSION}")
endforeach()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Install the library in the lib folder of the build directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

install(TARGETS ${xtargets}
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        #RUNTIME_DEPENDENCY_SET xdeps
)
#install(RUNTIME_DEPENDENCY_SET xdeps
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
#)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Recipe for installing the runtime dependencies.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#   The following dependency resolution approach yields the following CMake error:
#   ```bash
#   CMake Error at xrc/cmake_install.cmake:117 (file):
#     file Multiple conflicting paths found for libgcc_s.so.1:
#
#       /lib/x86_64-linux-gnu/libgcc_s.so.1
#       /usr/local/MATLAB/R2023b/bin/glnxa64/../../sys/os/glnxa64/libgcc_s.so.1
#   Call Stack (most recent call first):
#     cmake_install.cmake:52 (include)
#   ```
#   For now, the following dependencies checks are disabled until a better solution is found.
#   The major dependency that causes problems for MATLAB mex files is intel libirng.so.
#   But this dependency is now automatically linked statically at compile time using ``-static-intel`` flag.
#   The rest of the dependencies should be automatically copied to the install directory at the time of building libparamote.

# if (DLL_ENABLED AND DEPS_ENABLED)
#     install(CODE [[
#         function(install_library_with_deps LIBRARY)
#             file(INSTALL
#                 DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
#                 TYPE SHARED_LIBRARY
#                 FOLLOW_SYMLINK_CHAIN
#                 FILES "${LIBRARY}"
#             )
#             file(GET_RUNTIME_DEPENDENCIES
#                 LIBRARIES ${LIBRARY}
#                 RESOLVED_DEPENDENCIES_VAR resolved_deps
#                 UNRESOLVED_DEPENDENCIES_VAR unresolved_deps
#             )
#             foreach(FILE ${resolved_deps})
#                 if(NOT IS_SYMLINK ${FILE})
#                     install_library_with_deps(${FILE})
#                 endif()
#             endforeach()
#             foreach(FILE ${unresolved_deps})
#                 message(STATUS "Unresolved from ${LIBRARY}: ${FILE}")
#             endforeach()
#         endfunction()
#         file(GET_RUNTIME_DEPENDENCIES
#             LIBRARIES $<TARGET_SONAME_FILE:pm_sampling>
#             RESOLVED_DEPENDENCIES_VAR resolved_deps
#             UNRESOLVED_DEPENDENCIES_VAR unresolved_deps
#             POST_INCLUDE_REGEXES "libirng\\.so"
#             PRE_INCLUDE_REGEXES "libirng\\.so"
#             #POST_INCLUDE_REGEXES ".*[Ll][Ii][Bb][Ii][Rr][Nn][Gg].*"
#         )
#         foreach(FILE ${resolved_deps})
#             install_library_with_deps(${FILE})
#         endforeach()
#         foreach(FILE ${unresolved_deps})
#             message(STATUS "Unresolved: ${FILE}")
#         endforeach()
#     ]])
# endif()
