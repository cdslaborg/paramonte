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
#set(MEX_FLAGS "-v -nojvm")
#set(CFLAGS "-fPIC -shared")
#set(LDFLAGS "-fPIC -shared")

unset(rpath_prop)
if (APPLE)
    set(rpath_prop "@loader_path")
    set(LDFLAGS "${LDFLAGS}" "-static-libstdc++")
    set(LDFLAGS "${LDFLAGS}" "-Wl,-rpath,@rpath -Wl,-rpath,@loader_path")
elseif(UNIX)
    set(rpath_prop "$ORIGIN")
    set(LDFLAGS "${LDFLAGS}" "-static-libstdc++")
    set(LDFLAGS "${LDFLAGS}" "-Wl,-rpath,\$ORIGIN")
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
        target_compile_options("${xtarget}" BEFORE PUBLIC "${MEX_FLAGS}")
    endif()
    if (DEFINED LDFLAGS)
        target_link_options("${xtarget}" BEFORE PUBLIC "${LDFLAGS}")
    endif()
    set(xtargets "${xtargets}" "${xtarget}")
    add_dependencies(matlabmex "${xtarget}")
    get_target_property(MXLIB_INCLUDE_DIRECTORIES   "${xtarget}" INCLUDE_DIRECTORIES)
    get_target_property(MXLIB_COMPILE_DEFINITIONS   "${xtarget}" COMPILE_DEFINITIONS)
    get_target_property(MXLIB_COMPILE_OPTIONS       "${xtarget}" COMPILE_OPTIONS)
    get_target_property(MXLIB_COMPILE_FLAGS         "${xtarget}" COMPILE_FLAGS)
    get_target_property(MXLIB_LINK_DIRECTORIES  "${xtarget}" LINK_DIRECTORIES)
    get_target_property(MXLIB_INSTALL_NAME_DIR  "${xtarget}" INSTALL_NAME_DIR)
    get_target_property(MXLIB_LINK_OPTIONS      "${xtarget}" LINK_OPTIONS)
    get_target_property(MXLIB_LINK_FLAGS        "${xtarget}" LINK_FLAGS)
    message(NOTICE "${pmattn} MATLAB Mex target include directories : ${MXLIB_INCLUDE_DIRECTORIES}")
    message(NOTICE "${pmattn} MATLAB Mex target compile definitions : ${MXLIB_COMPILE_DEFINITIONS}")
    message(NOTICE "${pmattn} MATLAB Mex target compile options     : ${MXLIB_COMPILE_OPTIONS}")
    message(NOTICE "${pmattn} MATLAB Mex target compile flags       : ${MXLIB_COMPILE_FLAGS}")
    message(NOTICE "${pmattn} MATLAB Mex target link directories    : ${MXLIB_LINK_DIRECTORIES}")
    message(NOTICE "${pmattn} MATLAB Mex target install name dir    : ${MXLIB_INSTALL_NAME_DIR}")
    message(NOTICE "${pmattn} MATLAB Mex target link options        : ${MXLIB_LINK_OPTIONS}")
    message(NOTICE "${pmattn} MATLAB Mex target link flags          : ${MXLIB_LINK_FLAGS}")
    message(NOTICE "${pmattn} MATLAB Mex CMAKE INSTALL PATH         : ${paramonte_bld_lib_dir}")
endforeach()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Install the library in the lib folder of the build directory.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

install(TARGETS ${xtargets}
        CONFIGURATIONS "${CMAKE_BUILD_TYPE}"
        DESTINATION "${paramonte_bld_lib_dir}/"
        )